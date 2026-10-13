#include <string>
#include <cstdlib>
#include <stack>
#include <algorithm>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <cstring>
#include <unordered_set>
#include <array>
#include <sdsl/bit_vectors.hpp>
#include "include/Adaptive.hpp"
#include "include/LiaoParser.hpp"
#include "include/libsais64.h"


namespace {

class FmWaveletRank {
    static constexpr unsigned levels = 9;
    std::array<sdsl::bit_vector, levels> bits_;
    std::array<sdsl::rank_support_v<1>, levels> ranks_;
    std::array<uint32_t, levels> zeros_{};

public:
    void build(const std::vector<uint16_t>& values) {
        std::vector<uint16_t> current(values), next(values.size());
        for (unsigned level = 0; level < levels; ++level) {
            const unsigned shift = levels - 1 - level;
            bits_[level] = sdsl::bit_vector(values.size(), 0);
            uint32_t zero_count = 0;
            for (size_t i = 0; i < current.size(); ++i) {
                const bool bit = ((current[i] >> shift) & 1U) != 0;
                bits_[level][i] = bit;
                if (!bit) ++zero_count;
            }
            zeros_[level] = zero_count;
            size_t zero_position = 0, one_position = zero_count;
            for (uint16_t value : current) {
                if (((value >> shift) & 1U) == 0)
                    next[zero_position++] = value;
                else
                    next[one_position++] = value;
            }
            current.swap(next);
            ranks_[level] = sdsl::rank_support_v<1>(&bits_[level]);
        }
    }

    uint32_t rank(uint32_t position, uint16_t symbol) const {
        uint32_t left = 0, right = position;
        for (unsigned level = 0; level < levels; ++level) {
            const uint32_t left_ones = ranks_[level](left);
            const uint32_t right_ones = ranks_[level](right);
            const unsigned shift = levels - 1 - level;
            if (((symbol >> shift) & 1U) == 0) {
                left -= left_ones;
                right -= right_ones;
            } else {
                left = zeros_[level] + left_ones;
                right = zeros_[level] + right_ones;
            }
        }
        return right - left;
    }
};

} // namespace

BackwardEncoder::BackwardEncoder(WordBasedText* source)
    : n(source->Nchar), text(source) {
    codeStream = new uint32_t[n / 4 + 10];
    const std::string bytes = text->buffer.str();
    inputText = new unsigned char[n + 2]();
    std::memcpy(inputText, bytes.data(), n);
}

TextParser::TextParser(WordBasedText* w)
    : BackwardEncoder(w),
      langstr(new uchar[static_cast<size_t>(const_len) + 256 + depth]()),
      lang_capacity(static_cast<size_t>(const_len) + 256 + depth) {
    Lengths = new uint32_t[n];
    Offsets = new uint32_t[n];
    LangStrings[0].insert(lang);
}

TextParser::TextParser(WordBasedText* w, unsigned kmer_length,
                       size_t target_dictionary_size)
    : BackwardEncoder(w),
      langstr(new uchar[target_dictionary_size + depth]()),
      lang_capacity(target_dictionary_size + depth) {
    if (kmer_length == 0 || target_dictionary_size == 0 || n == 0)
        throw std::invalid_argument(
            "Uniform k-mer dictionary requires nonempty text, k, and size");

    bool present[256] = {false};
    for (uint32_t i = 0; i < n; ++i)
        present[inputText[i]] = true;

    size_t position = 0;
    for (unsigned value = 0; value < 256; ++value) {
        if (present[value]) {
            if (position == target_dictionary_size)
                throw std::invalid_argument(
                    "Target dictionary is smaller than the text alphabet");
            langstr[position++] = static_cast<uchar>(value);
        }
    }

    const size_t sampled_bytes = target_dictionary_size - position;
    const size_t sample_count =
        (sampled_bytes + kmer_length - 1) / kmer_length;
    const size_t last_start = n > kmer_length ? n - kmer_length : 0;

    for (size_t sample = 0; sample < sample_count; ++sample) {
        const size_t source = sample_count <= 1
            ? 0
            : static_cast<size_t>(
                  (static_cast<unsigned long long>(sample) * last_start) /
                  (sample_count - 1));
        for (unsigned q = 0;
             q < kmer_length && position < target_dictionary_size; ++q) {
            langstr[position++] = inputText[(source + q) % n];
        }
    }

    len = static_cast<int>(position);
    extra_len = 0;
    std::fill(std::begin(symbol_pos), std::end(symbol_pos), -1);
    for (int i = 0; i < len; ++i)
        if (symbol_pos[langstr[i]] == -1)
            symbol_pos[langstr[i]] = i;

    Lengths = new uint32_t[n];
    Offsets = new uint32_t[n];
}



TextParser::~TextParser() {
    delete[] codeStream;
    codeStream = nullptr;
    delete[] SA;
    delete[] LCP;
    delete RMQ;
    delete[] f;
    delete[] f1;
    delete[] langstr;
    delete[] Lengths;
    delete[] Offsets;
}

void TextParser::calc_imbalance() {
    for(auto const& i : stat[eulertigs_k-2]) {
        string vb = i.first.substr(0,i.first.length()-1);
        string ve = i.first.substr(1,i.first.length()-1);
        deBrujn.insert(make_pair(vb,arc{ve,0}));
        imbalance[vb]--;
        imbalance[ve]++;
    }
}

void TextParser::set_dictionary(const std::vector<uint8_t>& dictionary) {
    if (dictionary.empty() ||
        dictionary.size() > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Dictionary size is outside the supported range");
    if (dictionary.size() + depth > lang_capacity) {
        delete[] langstr;
        langstr = new uchar[dictionary.size() + depth]();
        lang_capacity = dictionary.size() + depth;
    } else {
        std::memset(langstr, 0, lang_capacity);
    }
    std::memcpy(langstr, dictionary.data(), dictionary.size());
    len = static_cast<int>(dictionary.size());
    extra_len = 0;
    blocked_dictionary_mode = false;
    common_dictionary_size = 0;
    regular_block_part_size = 0;
    dictionary_block_count = 1;
    text_block_size = 0;
    set_symbol_positions(langstr, dictionary.size());
}

void TextParser::build_eulertigs_dictionary(int order) {
    if (order < 3 || order > depth)
        throw std::invalid_argument("EulerTigs order is outside the supported range");
    eulertigs_k = order;
    for (auto& table : stat)
        table.clear();
    deBrujn.clear();
    imbalance.clear();
    cycle.clear();
    symbols.clear();

    // Only the vertices of order L-1 and edges of order L are required.
    std::array<bool, 256> text_symbols{};
    for (uint32_t i = 0; i < n; ++i) {
        text_symbols[inputText[i]] = true;
        if (i + order - 1 <= n)
            ++stat[order - 3][string(
                reinterpret_cast<char*>(inputText + i), order - 1)];
        if (i + order <= n)
            ++stat[order - 2][string(
                reinterpret_cast<char*>(inputText + i), order)];
    }
    for (unsigned value = 0; value < 256; ++value)
        if (text_symbols[value])
            symbols.insert(static_cast<unsigned char>(value));
    Nsymb = static_cast<int>(symbols.size());
    calc_imbalance();
    Eulerise();
    Euler();
    Spell();

    std::vector<uint8_t> dictionary(lang.begin(), lang.end());
    std::array<bool, 256> in_dictionary{};
    std::array<bool, 256> in_text{};
    for (uint8_t value : dictionary)
        in_dictionary[value] = true;
    for (uint32_t i = 0; i < n; ++i)
        in_text[inputText[i]] = true;
    for (unsigned value = 0; value < 256; ++value)
        if (in_text[value] && !in_dictionary[value])
            dictionary.push_back(static_cast<uint8_t>(value));
    set_dictionary(dictionary);
}

void TextParser::Eulerise() {
stack<string> st_minus,st_plus;
    for(auto const& i : stat[eulertigs_k-3]) {
        if(imbalance[i.first]>0) {
            if(!st_minus.empty()) {
                string s = st_minus.top();
                deBrujn.insert(make_pair(i.first,arc{s,1}));
                imbalance[i.first]--;
                imbalance[s]++;
                if(imbalance[s]==0)
                    st_minus.pop();
            }
            if(imbalance[i.first]>0)
                st_plus.push(i.first);
        } else if(imbalance[i.first]<0) {
            if(!st_plus.empty()) {
                string s = st_plus.top();
                deBrujn.insert(make_pair(s,arc{i.first,1}));
                imbalance[i.first]++;
                imbalance[s]--;
                if(imbalance[s]==0)
                    st_plus.pop();
            }
            if(imbalance[i.first]<0)
                st_minus.push(i.first);
        }
    }
    while(!st_minus.empty()) {
        string in = st_minus.top();
        string out = st_plus.top();
        deBrujn.insert(make_pair(out,arc{in,1}));
        imbalance[in]++;
        imbalance[out]--;
        if(imbalance[in]==0)
            st_minus.pop();
        if(imbalance[out]==0)
            st_plus.pop();
    }
}

void TextParser::Euler() {
string s((char*)inputText,eulertigs_k-1);
list<arc> cycle_;
    cycle.push_back(arc{s,0});
    auto first = cycle.begin();
    do {
        auto it = deBrujn.find(s);
        while(it!=deBrujn.end()) {
            cycle_.push_back(it->second);
            s = it->second.label;
            deBrujn.erase(it);
            it = deBrujn.find(s);
        }
        auto tmp = first;
        cycle.splice(++tmp, cycle_, cycle_.begin(), cycle_.end());
        while(first!=cycle.end() && deBrujn.find(first->label)==deBrujn.end())
            first++;
        if(first!=cycle.end())
            s = first->label;
    } while(first!=cycle.end());
}

void TextParser::Spell() {
auto it=cycle.begin();
    lang = it->label;
bool new_segm = false;
    for(it++;it!=cycle.end();it++) {
        if(!it->extra)
            if(new_segm) {
                lang += it->label;
                new_segm = false;
            } else
                lang += it->label.substr(it->label.length()-1,it->label.length()-1);
        else
            new_segm = true;
    }
}

void TextParser::build_suffix_array_dictionary(size_t dictionary_size,
                                               uint32_t max_suffix_length,
                                               int length_correction,
                                               const std::string& seed) {
    if (dictionary_size == 0 || n == 0 ||
        dictionary_size > static_cast<size_t>(std::numeric_limits<int>::max()))
        throw std::invalid_argument("Dictionary size is outside the supported range");
    if (max_suffix_length == 0)
        throw std::invalid_argument("Maximum suffix length must be positive");
    if (length_correction < 0)
        length_correction = n < 1000000 ? 0 : 1;
    if (length_correction != 0 && length_correction != 1)
        throw std::invalid_argument("Length correction b must be 0 or 1");

    bool alphabet_present[256] = {false};
    size_t alphabet_size = 0;
    symbols.clear();
    for (uint32_t position = 0; position < n; ++position) {
        const uchar symbol = inputText[position];
        symbols.insert(symbol);
        if (!alphabet_present[symbol]) {
            alphabet_present[symbol] = true;
            ++alphabet_size;
        }
    }
    if (seed.empty())
        throw std::invalid_argument("Dictionary seed must not be empty");
    if (dictionary_size < seed.size() + alphabet_size)
        throw std::invalid_argument(
            "Dictionary is too small for the seed and text alphabet");

    using DictionaryClock = std::chrono::high_resolution_clock;
    auto stage_begin = DictionaryClock::now();
    struct Interval { uint32_t left = 0, right = 0; };
    const uint32_t indexed_length = n + 1;
    std::vector<int64_t> reversed_text(indexed_length);
    for (uint32_t position = 0; position < n; ++position)
        reversed_text[position] =
            static_cast<int64_t>(inputText[n - 1 - position]) + 1;
    reversed_text[n] = 0; // unique sentinel

    std::vector<int64_t> text_sa(indexed_length);
    if (libsais64_long(reversed_text.data(), text_sa.data(), indexed_length,
                       257, 0) != 0)
        throw std::runtime_error("Reversed-text suffix-array construction failed");

    uint64_t counts[257] = {};
    for (int64_t symbol : reversed_text)
        ++counts[static_cast<unsigned>(symbol)];
    uint64_t cumulative[258] = {};
    for (unsigned symbol = 0; symbol < 257; ++symbol)
        cumulative[symbol + 1] = cumulative[symbol] + counts[symbol];

    std::vector<uint16_t> bwt(indexed_length);
    for (uint32_t rank = 0; rank < indexed_length; ++rank) {
        const uint32_t suffix = static_cast<uint32_t>(text_sa[rank]);
        const uint32_t preceding = suffix == 0 ? indexed_length - 1
                                                : suffix - 1;
        bwt[rank] = static_cast<uint16_t>(reversed_text[preceding]);
    }
    FmWaveletRank bwt_rank;
    bwt_rank.build(bwt);

    Interval one_byte[256];
    for (unsigned byte = 0; byte < 256; ++byte) {
        const unsigned symbol = byte + 1;
        one_byte[byte] = {
            static_cast<uint32_t>(cumulative[symbol]),
            static_cast<uint32_t>(cumulative[symbol + 1])};
    }
    reversed_text.clear();
    reversed_text.shrink_to_fit();
    text_sa.clear();
    text_sa.shrink_to_fit();
    bwt.clear();
    bwt.shrink_to_fit();
    auto stage_end = DictionaryClock::now();
    dictionary_stat_time_ms =
        std::chrono::duration<double, std::milli>(stage_end - stage_begin).count();

    struct ActiveSuffix {
        uint32_t length;
        Interval interval;
        uint64_t hash1;
        uint64_t hash2;
    };
    struct SubstringKey {
        uint64_t hash1, hash2;
        uint32_t length;
        bool operator==(const SubstringKey& other) const {
            return hash1 == other.hash1 && hash2 == other.hash2 &&
                   length == other.length;
        }
    };
    struct SubstringKeyHash {
        size_t operator()(const SubstringKey& key) const {
            uint64_t value = key.hash1 ^
                ((key.hash2 << 23) | (key.hash2 >> 41)) ^
                (static_cast<uint64_t>(key.length) * 0x9e3779b97f4a7c15ULL);
            value ^= value >> 30;
            value *= 0xbf58476d1ce4e5b9ULL;
            value ^= value >> 27;
            return static_cast<size_t>(value);
        }
    };
    constexpr uint64_t base1 = 0x9e3779b185ebca87ULL;
    constexpr uint64_t base2 = 0xc2b2ae3d27d4eb4fULL;

    auto refine = [&](const ActiveSuffix& suffix, uchar symbol) {
        const uint64_t mapped = static_cast<uint64_t>(symbol) + 1;
        return Interval{
            static_cast<uint32_t>(cumulative[mapped] +
                bwt_rank.rank(suffix.interval.left,
                              static_cast<uint16_t>(mapped))),
            static_cast<uint32_t>(cumulative[mapped] +
                bwt_rank.rank(suffix.interval.right,
                              static_cast<uint16_t>(mapped)))};
    };

    string dictionary;
    dictionary.reserve(dictionary_size);
    std::vector<ActiveSuffix> active;
    active.reserve(max_suffix_length);
    std::unordered_set<SubstringKey, SubstringKeyHash> represented;
    represented.reserve(dictionary_size * 8);

    // Add the new dictionary suffixes to the persistent represented set.
    auto append_symbol = [&](uchar symbol) {
        std::vector<ActiveSuffix> next;
        next.reserve(active.size() + 1);
        const Interval single = one_byte[symbol];
        if (single.left != single.right) {
            represented.insert({static_cast<uint64_t>(symbol) + 1,
                                static_cast<uint64_t>(symbol) + 1, 1});
            next.push_back({1, single, static_cast<uint64_t>(symbol) + 1,
                            static_cast<uint64_t>(symbol) + 1});
        }
        for (const ActiveSuffix& suffix : active) {
            if (suffix.length >= max_suffix_length)
                break;
            const Interval interval = refine(suffix, symbol);
            if (interval.left == interval.right)
                break;
            ActiveSuffix extended{
                suffix.length + 1, interval,
                suffix.hash1 * base1 + static_cast<uint64_t>(symbol) + 1,
                suffix.hash2 * base2 + static_cast<uint64_t>(symbol) + 1};
            represented.insert(
                {extended.hash1, extended.hash2, extended.length});
            next.push_back(extended);
        }
        dictionary.push_back(static_cast<char>(symbol));
        active.swap(next);
    };

    for (uchar symbol : seed)
        append_symbol(symbol);

    auto grow_to = [&](size_t target_size) {
        while (dictionary.size() < target_size) {
            uint64_t best_score = 0;
            uchar best_symbol = *symbols.begin();
            for (uchar candidate : symbols) {
                bool still_skipping = true;
                uint64_t score = 0;
                const Interval single = one_byte[candidate];
                const SubstringKey single_key{
                    static_cast<uint64_t>(candidate) + 1,
                    static_cast<uint64_t>(candidate) + 1, 1};
                if (single.left != single.right) {
                    if (represented.find(single_key) == represented.end()) {
                        still_skipping = false;
                        score += static_cast<uint64_t>(single.right - single.left) *
                                 static_cast<uint64_t>(1 - length_correction);
                    }
                }
                for (size_t p = 0; p < active.size(); ++p) {
                    const ActiveSuffix& suffix = active[p];
                    if (suffix.length >= max_suffix_length)
                        break;
                    const Interval interval = refine(suffix, candidate);
                    if (interval.left == interval.right)
                        break;
                    const SubstringKey key{
                        suffix.hash1 * base1 +
                            static_cast<uint64_t>(candidate) + 1,
                        suffix.hash2 * base2 +
                            static_cast<uint64_t>(candidate) + 1,
                        suffix.length + 1};
                    if (still_skipping && represented.find(key) != represented.end()) {
                        continue;
                    }
                    still_skipping = false;
                    score += static_cast<uint64_t>(interval.right - interval.left) *
                             static_cast<uint64_t>(suffix.length + 1 -
                                                   length_correction);
                }
                if (score > best_score) {
                    best_score = score;
                    best_symbol = candidate;
                }
            }
            if (best_score == 0) {
                size_t choice = static_cast<size_t>(rand()) % symbols.size();
                auto selected = symbols.begin();
                std::advance(selected, static_cast<long>(choice));
                best_symbol = *selected;
            }

            append_symbol(best_symbol);
        }
    };

    stage_begin = DictionaryClock::now();
    grow_to(dictionary_size - alphabet_size);
    for (unsigned value = 0; value < 256; ++value)
        if (alphabet_present[value] &&
            dictionary.find(static_cast<char>(value)) == string::npos)
            append_symbol(static_cast<uchar>(value));
    grow_to(dictionary_size);
    stage_end = DictionaryClock::now();
    dictionary_growth_time_ms =
        std::chrono::duration<double, std::milli>(stage_end - stage_begin).count();

    if (dictionary_size + depth > lang_capacity) {
        delete[] langstr;
        langstr = new uchar[dictionary_size + depth]();
        lang_capacity = dictionary_size + depth;
    } else {
        std::memset(langstr, 0, lang_capacity);
    }
    std::memcpy(langstr, dictionary.data(), dictionary_size);
    len = static_cast<int>(dictionary_size);
    extra_len = len - const_len;
    blocked_dictionary_mode = false;
    common_dictionary_size = 0;
    regular_block_part_size = 0;
    dictionary_block_count = 1;
    text_block_size = static_cast<int>(n);
    current_text_end = static_cast<int>(n);
    set_symbol_positions(langstr, dictionary_size);
}

void TextParser::set_symbol_positions(const uchar* dictionary, size_t size) {
    std::fill(std::begin(symbol_pos), std::end(symbol_pos), -1);
    for (size_t position = 0; position < size; ++position)
        if (symbol_pos[dictionary[position]] == -1)
            symbol_pos[dictionary[position]] = static_cast<int>(position);
}

void TextParser::process_text_liao(int repetitions, bool verbose) {
    if (repetitions <= 0)
        throw std::invalid_argument("RLZ parsing repetitions must be positive");
    if (blocked_dictionary_mode)
        throw std::invalid_argument(
            "The rlz-store parser does not support legacy split dictionaries");

    LiaoParser parser(static_cast<size_t>(n));
    auto begin = std::chrono::high_resolution_clock::now();
    parser.build_index(langstr, static_cast<size_t>(len));
    auto end = std::chrono::high_resolution_clock::now();
    liao_index_time_ms =
        std::chrono::duration<double, std::milli>(end - begin).count();

    begin = std::chrono::high_resolution_clock::now();
    for (int repetition = 0; repetition < repetitions; ++repetition)
        parser.parse(inputText, static_cast<size_t>(n));
    end = std::chrono::high_resolution_clock::now();
    liao_parse_time_ms =
        std::chrono::duration<double, std::milli>(end - begin).count() /
        repetitions;

    segn = static_cast<int>(parser.factor_count());
    std::copy_n(parser.offsets(), segn, Offsets);
    std::copy_n(parser.lengths(), segn, Lengths);
    max_len = 0;
    for (int i = 0; i < segn; ++i)
        max_len = std::max(max_len, static_cast<int>(Lengths[i]));

    if (verbose)
        cout << "=== Process text (rlz-store parser) ====\n"
             << "Segments:" << segn << " Code byte L=0 Avg length="
             << static_cast<double>(n) / segn << endl;
}

