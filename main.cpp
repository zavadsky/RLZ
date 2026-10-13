#include <algorithm>
#include <random>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <limits>
#include <filesystem>
#include <sstream>
#include <array>
#include <bit>
#include <cmath>
#include <numeric>
#include <unordered_map>
#include <unordered_set>
#include <stdint.h>
#include <zlib.h>
#include <chrono>
#include <dlfcn.h>
#include "include/WordBasedText.h"
#include "include/Adaptive.hpp"
#include "include/Rank.hpp"
#include "include/RLZ_Rank_sdsl.hpp"
#include "include/LiaoParser.hpp"

namespace {

using BenchmarkClock = std::chrono::high_resolution_clock;

double elapsed_ms(BenchmarkClock::time_point begin,
                  BenchmarkClock::time_point end) {
    return std::chrono::duration<double, std::milli>(end - begin).count();
}

// Prevent the optimizer from removing timed extraction loops.
volatile uint64_t benchmark_sink = 0;

struct SizeBenchmarkInput {
    std::string file;
    std::string seed = "the ";
    size_t segment_length = 0;
    uint32_t max_suffix_length = 0;
};

// Read the common experiment description used by all three output modes.
std::vector<SizeBenchmarkInput> read_size_benchmark_inputs(
    const std::string& csv_name) {
    std::ifstream input(csv_name);
    if (!input)
        throw std::runtime_error("Cannot open file list: " + csv_name);
    const std::filesystem::path base =
        std::filesystem::absolute(csv_name).parent_path();
    std::vector<SizeBenchmarkInput> result;
    std::string line;
    size_t line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (!line.empty() && line.back() == '\r')
            line.pop_back();
        if (line.empty())
            continue;
        std::array<std::string, 4> fields{"", "the ", "300", "24"};
        size_t begin = 0;
        for (size_t field = 0; field < fields.size() && begin <= line.size(); ++field) {
            const size_t end = line.find_first_of(";,\t", begin);
            fields[field] = line.substr(begin, end - begin);
            if (end == std::string::npos)
                break;
            begin = end + 1;
        }
        if (fields[0] == "file" || fields[0] == "File")
            continue;
        std::filesystem::path path(fields[0]);
        if (path.is_relative())
            path = base / path;
        SizeBenchmarkInput row;
        row.file = path.lexically_normal().string();
        row.seed = fields[1];
        if (row.seed.empty())
            throw std::runtime_error("Seed must not be empty at line " +
                                     std::to_string(line_number));
        try {
            row.segment_length = static_cast<size_t>(std::stoull(fields[2]));
            row.max_suffix_length =
                static_cast<uint32_t>(std::stoul(fields[3]));
        } catch (const std::exception&) {
            throw std::runtime_error("Invalid numeric field at line " +
                                     std::to_string(line_number));
        }
        if (row.segment_length < 16)
            throw std::runtime_error("Segment length must be at least 16");
        if (row.max_suffix_length == 0)
            throw std::runtime_error("L must be positive");
        result.push_back(std::move(row));
    }
    if (result.empty())
        throw std::runtime_error("The file list is empty");
    return result;
}

constexpr size_t liao_k = 16;
constexpr uint32_t liao_seed = 2;
using LiaoKmer = std::array<uint8_t, liao_k>;

struct LiaoKmerHash {
    size_t operator()(const LiaoKmer& word) const noexcept {
        uint64_t hash = 14695981039346656037ULL;
        for (uint8_t value : word) {
            hash ^= value;
            hash *= 1099511628211ULL;
        }
        return static_cast<size_t>(hash);
    }
};

LiaoKmer liao_word_at(const uint8_t* text, size_t position) {
    LiaoKmer word;
    std::copy_n(text + position, liao_k, word.data());
    return word;
}

std::vector<uint8_t> build_liao_dictionary(const uint8_t* data, size_t n,
                                           size_t budget,
                                           size_t segment,
                                           double scoring_power) {
    if (n < liao_k || segment < liao_k || budget == 0)
        throw std::invalid_argument("Text, segment, or dictionary is too small for Liao LMC");
    const size_t full_segments = n / segment;
    const size_t selected_count = (budget + segment - 1) / segment;
    if (full_segments == 0 || selected_count > full_segments)
        throw std::invalid_argument("Not enough full text segments for Liao LMC");

    std::array<bool, 256> input_symbols{};
    for (size_t i = 0; i < n; ++i)
        input_symbols[data[i]] = true;
    const size_t alphabet_size = static_cast<size_t>(
        std::count(input_symbols.begin(), input_symbols.end(), true));
    if (budget < alphabet_size)
        throw std::invalid_argument("Dictionary budget is smaller than the text alphabet");

    const size_t threshold = std::max<size_t>(
        1, std::min<size_t>(256, (n / budget) / 2));
    const size_t population = n - liao_k + 1;
    const size_t reservoir_size = (population + threshold - 1) / threshold;
    std::vector<LiaoKmer> reservoir;
    reservoir.reserve(reservoir_size);
    std::mt19937_64 sampler(liao_seed);
    for (size_t position = 0; position < population; ++position) {
        const LiaoKmer word = liao_word_at(data, position);
        if (position < reservoir_size) {
            reservoir.push_back(word);
        } else {
            const size_t slot =
                std::uniform_int_distribution<size_t>(0, position)(sampler);
            if (slot < reservoir_size)
                reservoir[slot] = word;
        }
    }
    std::unordered_map<LiaoKmer, uint32_t, LiaoKmerHash> frequencies;
    frequencies.reserve(reservoir.size());
    for (const LiaoKmer& word : reservoir)
        ++frequencies[word];

    std::vector<size_t> epochs(selected_count);
    std::iota(epochs.begin(), epochs.end(), 0);
    std::mt19937_64 order(liao_seed);
    std::shuffle(epochs.begin(), epochs.end(), order);
    std::unordered_set<LiaoKmer, LiaoKmerHash> covered;
    std::vector<size_t> chosen;
    chosen.reserve(selected_count);
    for (size_t epoch : epochs) {
        const size_t begin = epoch * full_segments / selected_count;
        const size_t end = (epoch + 1) * full_segments / selected_count;
        size_t best = begin;
        double best_score = -1.0;
        for (size_t candidate = begin; candidate < end; ++candidate) {
            std::unordered_set<LiaoKmer, LiaoKmerHash> local;
            double score = 0.0;
            const size_t start = candidate * segment;
            for (size_t offset = 0; offset + liao_k <= segment; ++offset) {
                const LiaoKmer word = liao_word_at(data, start + offset);
                if (covered.find(word) != covered.end() ||
                    !local.insert(word).second)
                    continue;
                const auto found = frequencies.find(word);
                if (found != frequencies.end()) {
                    if (scoring_power == 0.0)
                        score += 1.0;
                    else if (scoring_power == 0.5)
                        score += std::sqrt(static_cast<double>(found->second));
                    else
                        score += static_cast<double>(found->second);
                }
            }
            if (score > best_score) {
                best_score = score;
                best = candidate;
            }
        }
        chosen.push_back(best);
        const size_t start = best * segment;
        for (size_t offset = 0; offset + liao_k <= segment; ++offset)
            covered.insert(liao_word_at(data, start + offset));
    }
    std::sort(chosen.begin(), chosen.end());

    std::vector<uint8_t> dictionary;
    dictionary.reserve(budget);
    size_t remaining = budget;
    for (size_t index : chosen) {
        const size_t count = std::min(segment, remaining);
        dictionary.insert(dictionary.end(), data + index * segment,
                          data + index * segment + count);
        remaining -= count;
    }
    std::array<bool, 256> retained{};
    for (size_t i = 0; i < budget - alphabet_size; ++i)
        retained[dictionary[i]] = true;
    std::vector<uint8_t> missing;
    for (size_t value = 0; value < 256; ++value)
        if (input_symbols[value] && !retained[value])
            missing.push_back(static_cast<uint8_t>(value));
    dictionary.resize(budget - missing.size());
    dictionary.insert(dictionary.end(), missing.begin(), missing.end());
    return dictionary;
}

// Return the size of one independently Zlib-compressed byte stream.
uint64_t zlib_encoded_size(const void* source, size_t bytes) {
    if (bytes == 0)
        return 0;
    z_stream stream{};
    if (deflateInit2(&stream, 9, Z_DEFLATED, 15, 9,
                     Z_DEFAULT_STRATEGY) != Z_OK)
        throw std::runtime_error("deflateInit2 failed");
    std::vector<uint8_t> output(deflateBound(&stream, bytes));
    stream.next_in = reinterpret_cast<Bytef*>(const_cast<void*>(source));
    stream.avail_in = static_cast<uInt>(bytes);
    stream.next_out = output.data();
    stream.avail_out = static_cast<uInt>(output.size());
    const int status = deflate(&stream, Z_FINISH);
    const uint64_t result = stream.total_out;
    deflateEnd(&stream);
    if (status != Z_STREAM_END)
        throw std::runtime_error("Zlib compression failed");
    return result;
}

unsigned packed_bits_needed(uint64_t maximum) {
    return maximum == 0 ? 1U : 64U - std::countl_zero(maximum);
}

uint64_t sdsl_packed_vector_size(size_t count, uint64_t maximum) {
    const uint64_t bits = static_cast<uint64_t>(count) *
                          packed_bits_needed(maximum);
    return 9 + ((bits + 63) / 64) * 8;
}

uint64_t liao_archive_size(const uint8_t* text, size_t text_size,
                           LiaoParser& parser, size_t block_size) {
    constexpr uint32_t literal_threshold = 3;
    const size_t block_count = (text_size + block_size - 1) / block_size;
    uint64_t streams = 0;
    uint64_t stream_position_bits = 0;
    uint64_t maximum_stream_position_bits = 0;
    uint64_t maximum_factors = 0;

    for (size_t block_id = 0; block_id < block_count; ++block_id) {
        const size_t begin = block_id * block_size;
        const size_t size = std::min(block_size, text_size - begin);
        const size_t greedy_count = parser.parse(text + begin, size, false);
        const uint32_t* greedy_offsets = parser.offsets();
        const uint32_t* greedy_lengths = parser.lengths();
        std::vector<uint32_t> lengths;
        std::vector<uint32_t> offsets;
        std::vector<uint8_t> literals;
        lengths.reserve(greedy_count);
        offsets.reserve(greedy_count);
        size_t text_position = 0;
        bool previous_literal = false;
        for (size_t i = 0; i < greedy_count; ++i) {
            const uint32_t length = greedy_lengths[i];
            if (length <= literal_threshold) {
                if (previous_literal &&
                    lengths.back() + length <= literal_threshold)
                    lengths.back() += length;
                else
                    lengths.push_back(length);
                literals.insert(literals.end(), text + begin + text_position,
                                text + begin + text_position + length);
                previous_literal = true;
            } else {
                lengths.push_back(length);
                offsets.push_back(greedy_offsets[i]);
                previous_literal = false;
            }
            text_position += length;
        }
        for (uint32_t& length : lengths)
            --length;
        uint64_t block_bytes = 4 + zlib_encoded_size(
            lengths.data(), lengths.size() * sizeof(uint32_t));
        if (!literals.empty())
            block_bytes += 4 + zlib_encoded_size(literals.data(),
                                                  literals.size());
        if (!offsets.empty())
            block_bytes += 4 + zlib_encoded_size(
                offsets.data(), offsets.size() * sizeof(uint32_t));
        maximum_stream_position_bits = stream_position_bits;
        stream_position_bits += 8 * block_bytes;
        streams += block_bytes;
        maximum_factors = std::max<uint64_t>(maximum_factors,
                                             lengths.size());
    }
    return streams +
        sdsl_packed_vector_size(block_count, maximum_stream_position_bits) +
        sdsl_packed_vector_size(block_count, maximum_factors);
}

// Quote a field only when required by RFC-style CSV syntax.
std::string csv_field(const std::string& value) {
    if (value.find_first_of(",\"\r\n") == std::string::npos)
        return value;
    std::string result = "\"";
    for (char c : value) {
        result.push_back(c);
        if (c == '"')
            result.push_back('"');
    }
    result.push_back('"');
    return result;
}

struct QualityResult {
    size_t dictionary_size = 0;
    int segments = 0;
    double average_segment_length = 0.0;
    double dictionary_ms = 0.0;
};

QualityResult evaluate_dictionary_quality(TextParser& parser,
                                          const WordBasedText& text,
                                          double dictionary_ms) {
    parser.process_text_liao(1, false);
    if (parser.segments_number() <= 0)
        throw std::runtime_error("RLZ factorization produced no phrases");
    return {static_cast<size_t>(parser.dictionary_size()),
            parser.segments_number(),
            static_cast<double>(text.Nchar) / parser.segments_number(),
            dictionary_ms};
}


int run_dictionary_quality_benchmark(const std::string& list_name) {
    const auto inputs = read_size_benchmark_inputs(list_name);
    const std::filesystem::path output_name =
        std::filesystem::absolute(list_name).parent_path() /
        "table1.csv";
    std::ofstream csv(output_name);
    if (!csv)
        throw std::runtime_error("Cannot create " + output_name.string());
    csv << "File,Dictionary size,"
           "EulerTigs average phrase length,EulerTigs construction time (ms),"
           "Uniform average phrase length,Uniform construction time (ms),"
           "LMC average phrase length,LMC construction time (ms),LMC ell,"
           "Our method average phrase length,Our method construction time (ms),EulerTigs k\n";
    csv << std::setprecision(12);

    struct MethodResult {
        bool valid = false;
        double average = 0.0;
        double construction_ms = 0.0;
    };

    auto show_result = [](const char* method, const QualityResult& result) {
        std::cout << "  " << method << ": dictionary="
                  << result.dictionary_size << " segments=" << result.segments
                  << " average length=" << result.average_segment_length
                  << " construction time=" << result.dictionary_ms << " ms"
                  << std::endl;
    };

    for (const SizeBenchmarkInput& input : inputs) {
        const std::string& file = input.file;
        std::cout << "Dictionary-quality benchmark: " << file << std::endl;
        WordBasedText text(file, 1);
        const size_t minimum_eulertigs_size = static_cast<size_t>(
            std::floor(0.014L * static_cast<long double>(text.Nchar))) + 1;
        int eulertigs_order = 3;
        size_t target = 0;
        for (;; ++eulertigs_order) {
            TextParser probe(&text);
            probe.build_eulertigs_dictionary(eulertigs_order);
            target = static_cast<size_t>(probe.dictionary_size());
            if (target >= minimum_eulertigs_size)
                break;
            if (eulertigs_order >= 64)
                throw std::runtime_error(
                    "No EulerTigs dictionary reached 1.4% of the text");
        }
        constexpr size_t block_size = 300;
        std::cout << "  text=" << text.Nchar << " target=" << target
                  << " EulerTigs k=" << eulertigs_order
                  << " block=" << block_size << " L=24" << std::endl;

        MethodResult eulertigs, uniform, ours;
        MethodResult liao_l0, liao_l05, liao_l1;
        try {
            const auto begin = BenchmarkClock::now();
            TextParser parser(&text);
            parser.build_eulertigs_dictionary(eulertigs_order);
            const double dictionary_ms = elapsed_ms(begin, BenchmarkClock::now());
            const QualityResult result =
                evaluate_dictionary_quality(parser, text, dictionary_ms);
            eulertigs = {true, result.average_segment_length,
                         result.dictionary_ms};
            const std::string label =
                "EulerTigs (L=" + std::to_string(eulertigs_order) + ")";
            show_result(label.c_str(), result);
        } catch (const std::exception& error) {
            std::cerr << "  EulerTigs failed: " << error.what() << std::endl;
        }

        try {
            const auto begin = BenchmarkClock::now();
            TextParser parser(&text, block_size, target);
            const double dictionary_ms = elapsed_ms(begin, BenchmarkClock::now());
            const QualityResult result =
                evaluate_dictionary_quality(parser, text, dictionary_ms);
            uniform = {true, result.average_segment_length,
                       result.dictionary_ms};
            show_result("Uniform sweep", result);
        } catch (const std::exception& error) {
            std::cerr << "  Uniform sweep failed: " << error.what() << std::endl;
        }

        const auto run_liao = [&](size_t segment, double scoring_power,
                                  const char* label) -> MethodResult {
            try {
                TextParser parser(&text);
                const auto begin = BenchmarkClock::now();
                const std::vector<uint8_t> dictionary = build_liao_dictionary(
                    parser.text_data(), static_cast<size_t>(text.Nchar), target,
                    segment, scoring_power);
                parser.set_dictionary(dictionary);
                const double dictionary_ms =
                    elapsed_ms(begin, BenchmarkClock::now());
                const QualityResult result =
                    evaluate_dictionary_quality(parser, text, dictionary_ms);
                show_result(label, result);
                return {true, result.average_segment_length,
                        result.dictionary_ms};
            } catch (const std::exception& error) {
                std::cerr << "  " << label << " failed: " << error.what()
                          << std::endl;
                return {};
            }
        };

        const std::string liao_l0_label =
            "Liao l0 s=" + std::to_string(block_size);
        const std::string liao_l05_label =
            "Liao l0.5 s=" + std::to_string(block_size);
        const std::string liao_l1_label =
            "Liao l1 s=" + std::to_string(block_size);
        liao_l0 = run_liao(block_size, 0.0, liao_l0_label.c_str());
        liao_l05 = run_liao(block_size, 0.5, liao_l05_label.c_str());
        liao_l1 = run_liao(block_size, 1.0, liao_l1_label.c_str());

        try {
            TextParser parser(&text);
            const auto begin = BenchmarkClock::now();
            parser.build_suffix_array_dictionary(target, 24, -1, input.seed);
            const double dictionary_ms = elapsed_ms(begin, BenchmarkClock::now());
            const QualityResult result =
                evaluate_dictionary_quality(parser, text, dictionary_ms);
            ours = {true, result.average_segment_length,
                    result.dictionary_ms};
            show_result("Our v62", result);
        } catch (const std::exception& error) {
            std::cerr << "  Our v62 failed: " << error.what() << std::endl;
        }

        const MethodResult* best_liao = &liao_l0;
        double best_liao_ell = 0.0;
        if (liao_l05.valid && (!best_liao->valid ||
                              liao_l05.average > best_liao->average)) {
            best_liao = &liao_l05;
            best_liao_ell = 0.5;
        }
        if (liao_l1.valid && (!best_liao->valid ||
                             liao_l1.average > best_liao->average)) {
            best_liao = &liao_l1;
            best_liao_ell = 1.0;
        }

        csv << csv_field(std::filesystem::path(file).filename().string())
            << ',' << target;
        const auto write_method = [&](const MethodResult& result) {
            if (result.valid)
                csv << ',' << result.average << ',' << result.construction_ms;
            else
                csv << ",,";
        };
        write_method(eulertigs);
        write_method(uniform);
        write_method(*best_liao);
        csv << ',' << best_liao_ell;
        write_method(ours);
        csv << ',' << eulertigs_order << '\n';
    }
    std::cout << "Dictionary-quality report: " << output_name.string()
              << std::endl;
    return 0;
}

std::array<uint64_t, 4> four_representation_sizes(TextParser& parser,
                                                  WordBasedText& text) {
    parser.process_text_liao(1, false);
    RLZ_Rank pasta_fixed(
        parser.lengths_data(), parser.offsets_data(), parser.dictionary_data(),
        parser.segments_number(), text.Nchar, parser.dictionary_size(), false,
        0, 0, 1, false);
    RLZ_Rank pasta_vbyte(
        parser.lengths_data(), parser.offsets_data(), parser.dictionary_data(),
        parser.segments_number(), text.Nchar, parser.dictionary_size(), false,
        0, 0, 1, true);
    RLZ_Rank_ef ef_fixed(
        parser.lengths_data(), parser.offsets_data(), parser.dictionary_data(),
        parser.segments_number(), text.Nchar, parser.dictionary_size(),
        0, 0, 1, false);
    RLZ_Rank_ef ef_vbyte(
        parser.lengths_data(), parser.offsets_data(), parser.dictionary_data(),
        parser.segments_number(), text.Nchar, parser.dictionary_size(),
        0, 0, 1, true);
    const uint64_t dictionary = parser.dictionary_size();
    return {
        dictionary + pasta_fixed.space_usage() +
            pasta_fixed.packed_offsets_size(),
        dictionary + pasta_vbyte.space_usage() +
            pasta_vbyte.packed_offsets_size(),
        dictionary + ef_fixed.space_usage() + ef_fixed.packed_offsets_size(),
        dictionary + ef_vbyte.space_usage() + ef_vbyte.packed_offsets_size()
    };
}

int run_structure_size_benchmark(const std::string& list_name) {
    const auto inputs = read_size_benchmark_inputs(list_name);
    const std::filesystem::path output_name =
        std::filesystem::absolute(list_name).parent_path() /
        "table2.csv";
    std::ofstream csv(output_name);
    if (!csv)
        throw std::runtime_error("Cannot create " + output_name.string());
    csv << "file,dictionary size,PASTA + fixed offsets,"
           "PASTA + VByte offsets,Elias-Fano + fixed offsets,"
           "Elias-Fano + VByte offsets,Uniform + PASTA + fixed offsets,"
           "Uniform + Elias-Fano + VByte offsets,Liao complete Zlib scheme\n";

    for (const SizeBenchmarkInput& input : inputs) {
        WordBasedText text(input.file, 1);
        const size_t text_size = static_cast<size_t>(text.Nchar);
        size_t dictionary_size = 1;
        while (static_cast<long double>(dictionary_size) <=
               0.014L * text_size) {
            if (dictionary_size >
                std::numeric_limits<size_t>::max() / 2)
                throw std::overflow_error("Dictionary size overflow");
            dictionary_size *= 2;
        }
        std::cout << "Structure-size benchmark: " << input.file
                  << " text=" << text_size
                  << " dictionary=" << dictionary_size
                  << " segment=" << input.segment_length
                  << " L=" << input.max_suffix_length << std::endl;

        std::array<uint64_t, 4> our_sizes;
        {
            TextParser parser(&text);
            parser.build_suffix_array_dictionary(dictionary_size,
                                                 input.max_suffix_length, -1,
                                                 input.seed);
            our_sizes = four_representation_sizes(parser, text);
        }

        uint64_t uniform_size = 0;
        uint64_t uniform_ef_vbyte_size = 0;
        {
            TextParser parser(&text, input.segment_length, dictionary_size);
            const auto uniform_sizes = four_representation_sizes(parser, text);
            uniform_size = uniform_sizes[0];
            uniform_ef_vbyte_size = uniform_sizes[3];
        }

        uint64_t liao_size = std::numeric_limits<uint64_t>::max();
        size_t best_liao_block = 0;
        {
            TextParser builder(&text);
            const std::vector<uint8_t> dictionary = build_liao_dictionary(
                builder.text_data(), text_size, dictionary_size,
                input.segment_length, 0.5);
            LiaoParser parser(text_size);
            // The complete Liao baseline uses its original factor selector.
            parser.build_index(dictionary.data(), dictionary.size(), false);
            for (size_t block_size = 100;; block_size *= 2) {
                const uint64_t candidate = dictionary_size +
                    liao_archive_size(builder.text_data(), text_size,
                                      parser, block_size);
                if (candidate < liao_size) {
                    liao_size = candidate;
                    best_liao_block = block_size;
                }
                if (block_size >= text_size)
                    break;
                if (block_size > std::numeric_limits<size_t>::max() / 2)
                    throw std::overflow_error("Liao block size overflow");
            }
        }

        std::cout << "  PASTA/fixed=" << our_sizes[0]
                  << " PASTA/VByte=" << our_sizes[1]
                  << " EF/fixed=" << our_sizes[2]
                  << " EF/VByte=" << our_sizes[3]
                  << " uniform=" << uniform_size
                  << " uniform EF/VByte=" << uniform_ef_vbyte_size
                  << " Liao=" << liao_size
                  << " (best block=" << best_liao_block << ")" << std::endl;
        csv << csv_field(std::filesystem::path(input.file)
                             .filename().string())
            << ',' << dictionary_size;
        for (uint64_t size : our_sizes)
            csv << ',' << size;
        csv << ',' << uniform_size << ',' << uniform_ef_vbyte_size
            << ',' << liao_size << '\n';
        csv.flush();
    }
    std::cout << "Structure-size report: " << output_name.string()
              << std::endl;
    return 0;
}

struct ZstdApi {
    using CompressBound = size_t (*)(size_t);
    using Compress = size_t (*)(void*, size_t, const void*, size_t, int);
    using Decompress = size_t (*)(void*, size_t, const void*, size_t);
    using IsError = unsigned (*)(size_t);
    using ErrorName = const char* (*)(size_t);
    void* handle = nullptr;
    CompressBound compress_bound = nullptr;
    Compress compress = nullptr;
    Decompress decompress = nullptr;
    IsError is_error = nullptr;
    ErrorName error_name = nullptr;

    ZstdApi() {
        handle = dlopen("libzstd.so.1", RTLD_NOW | RTLD_LOCAL);
        if (!handle)
            handle = dlopen("libzstd.so", RTLD_NOW | RTLD_LOCAL);
        if (!handle)
            throw std::runtime_error("libzstd is required for Figure 1 mode");
        compress_bound = reinterpret_cast<CompressBound>(dlsym(handle, "ZSTD_compressBound"));
        compress = reinterpret_cast<Compress>(dlsym(handle, "ZSTD_compress"));
        decompress = reinterpret_cast<Decompress>(dlsym(handle, "ZSTD_decompress"));
        is_error = reinterpret_cast<IsError>(dlsym(handle, "ZSTD_isError"));
        error_name = reinterpret_cast<ErrorName>(dlsym(handle, "ZSTD_getErrorName"));
        if (!compress_bound || !compress || !decompress || !is_error || !error_name)
            throw std::runtime_error("The installed libzstd lacks required symbols");
    }
    ~ZstdApi() { if (handle) dlclose(handle); }
    void check(size_t code, const char* operation) const {
        if (is_error(code))
            throw std::runtime_error(std::string(operation) + ": " + error_name(code));
    }
};

struct ZstdBlockArchive {
    size_t block_size = 0;
    size_t text_size = 0;
    uint64_t format_overhead = 0;
    std::vector<std::vector<uint8_t>> frames;

    uint64_t size() const {
        uint64_t result = format_overhead;
        for (const auto& frame : frames)
            result += frame.size();
        return result;
    }
};

ZstdBlockArchive build_zstd_blocks(const ZstdApi& zstd, const uint8_t* text,
                                   size_t n, size_t block_size,
                                   bool seekable_format) {
    ZstdBlockArchive result;
    result.block_size = std::max<size_t>(1, block_size);
    result.text_size = n;
    const size_t blocks = (n + result.block_size - 1) / result.block_size;
    // Seekable Zstd uses an 8-byte skippable-frame header, an 8-byte entry
    // per frame without checksums, and a 9-byte footer.
    result.format_overhead = seekable_format ? 17 + 8 * blocks : 96 + 8 * blocks;
    result.frames.reserve(blocks);
    for (size_t block = 0; block < blocks; ++block) {
        const size_t begin = block * result.block_size;
        const size_t bytes = std::min(result.block_size, n - begin);
        std::vector<uint8_t> compressed(zstd.compress_bound(bytes));
        const size_t written = zstd.compress(compressed.data(), compressed.size(),
                                             text + begin, bytes, 22);
        zstd.check(written, "ZSTD_compress");
        compressed.resize(written);
        result.frames.push_back(std::move(compressed));
    }
    return result;
}

ZstdBlockArchive closest_seekable_archive(const ZstdApi& zstd,
                                           const uint8_t* text, size_t n,
                                           uint64_t target) {
    const uint64_t tolerance = (target + 499) / 200; // 0.5%
    ZstdBlockArchive smallest_frames =
        build_zstd_blocks(zstd, text, n, 64, true);
    ZstdBlockArchive single_frame =
        build_zstd_blocks(zstd, text, n, std::max<size_t>(1, n), true);
    const uint64_t largest_size = smallest_frames.size();
    const uint64_t smallest_size = single_frame.size();
    if (target + tolerance < smallest_size) {
        throw std::runtime_error(
            "Target cannot be reached: even a single Zstd frame is " +
            std::to_string(100.0 * static_cast<double>(smallest_size - target) /
                           target) +
            "% larger than the target");
    }
    if (target > largest_size + tolerance) {
        throw std::runtime_error(
            "Target cannot be reached: 64-byte Zstd frames are still " +
            std::to_string(100.0 * static_cast<double>(target - largest_size) /
                           target) +
            "% smaller than the target");
    }

    size_t low = 64, high = std::max<size_t>(64, n);
    ZstdBlockArchive best = std::move(single_frame);
    uint64_t best_distance = best.size() > target
        ? best.size() - target : target - best.size();
    const uint64_t small_frame_distance = largest_size > target
        ? largest_size - target : target - largest_size;
    if (small_frame_distance < best_distance) {
        best_distance = small_frame_distance;
        best = std::move(smallest_frames);
    }
    if (best_distance <= tolerance)
        return best;
    auto consider = [&](size_t block) {
        ZstdBlockArchive candidate = build_zstd_blocks(zstd, text, n, block, true);
        const uint64_t distance = candidate.size() > target
            ? candidate.size() - target : target - candidate.size();
        if (distance < best_distance) {
            best_distance = distance;
            best = std::move(candidate);
        }
    };
    // Compression normally improves with frame size, enabling binary search.
    while (low <= high) {
        const size_t mid = low + (high - low) / 2;
        ZstdBlockArchive candidate = build_zstd_blocks(zstd, text, n, mid, true);
        const uint64_t candidate_size = candidate.size();
        const uint64_t distance = candidate_size > target
            ? candidate_size - target : target - candidate_size;
        if (distance < best_distance) {
            best_distance = distance;
            best = std::move(candidate);
        }
        if (best_distance <= tolerance)
            return best;
        // Larger frames normally compress better and therefore reduce size.
        if (candidate_size > target)
            low = mid + 1;
        else {
            if (mid == 0) break;
            high = mid - 1;
        }
    }
    const size_t centre = best.block_size;
    const size_t radius = std::max<size_t>(16, centre / 200);
    const size_t step = std::max<size_t>(1, radius / 16);
    for (size_t block = centre > radius ? centre - radius : 1;
         block <= std::min(n, centre + radius); block += step) {
        consider(block);
        if (best_distance <= tolerance)
            return best;
    }
    throw std::runtime_error(
        "The target is attainable, but the frame-size search did not reach "
        "0.2% precision; closest tested relative error is " +
        std::to_string(100.0 * static_cast<double>(best_distance) / target) +
        "%");
}

ZstdBlockArchive closest_squashfs_archive(const ZstdApi& zstd,
                                          const uint8_t* text, size_t n,
                                          uint64_t target,
                                          size_t seekable_frame_size) {
    ZstdBlockArchive best;
    uint64_t best_distance = std::numeric_limits<uint64_t>::max();
    constexpr size_t minimum_block = 4096;
    constexpr size_t maximum_block = 1024 * 1024;
    size_t lower = minimum_block;
    while (lower < maximum_block && lower * 2 <= seekable_frame_size)
        lower *= 2;
    const size_t upper = std::min(maximum_block, lower * 2);

    const auto consider = [&](size_t block) {
        ZstdBlockArchive candidate = build_zstd_blocks(zstd, text, n, block, false);
        const uint64_t distance = candidate.size() > target
            ? candidate.size() - target : target - candidate.size();
        if (distance < best_distance) {
            best_distance = distance;
            best = std::move(candidate);
        }
    };
    consider(lower);
    if (upper != lower)
        consider(upper);
    return best;
}

std::vector<double> benchmark_block_archive(const ZstdApi& zstd,
                                            const ZstdBlockArchive& archive,
                                            const uint8_t* original) {
    constexpr int queries = 10000;
    std::vector<double> times;
    std::mt19937 rng(12345);
    std::vector<size_t> positions(queries);
    std::vector<uint8_t> decoded(archive.block_size);
    std::vector<uint8_t> output;
    for (size_t length = 8; length <= 4096; length *= 2) {
        if (length > archive.text_size) { times.push_back(0.0); continue; }
        output.resize(length);
        std::uniform_int_distribution<size_t> distribution(0, archive.text_size - length);
        for (size_t q = 0; q < queries; ++q)
            positions[q] = distribution(rng);
        auto extract_query = [&](size_t query_begin) {
            const size_t query_end = query_begin + length;
            size_t out = 0;
            for (size_t block = query_begin / archive.block_size;
                 block <= (query_end - 1) / archive.block_size; ++block) {
                const size_t plain_begin = block * archive.block_size;
                const size_t plain_size = std::min(archive.block_size,
                                                   archive.text_size - plain_begin);
                const auto& frame = archive.frames[block];
                const size_t produced = zstd.decompress(decoded.data(), decoded.size(),
                                                        frame.data(), frame.size());
                zstd.check(produced, "ZSTD_decompress");
                if (produced != plain_size)
                    throw std::runtime_error("Unexpected decompressed block size");
                const size_t from = std::max(query_begin, plain_begin) - plain_begin;
                const size_t to = std::min(query_end, plain_begin + plain_size) - plain_begin;
                std::copy(decoded.begin() + from, decoded.begin() + to,
                          output.begin() + out);
                out += to - from;
            }
        };
        uint64_t checksum = 0;
        const auto begin_time = BenchmarkClock::now();
        for (size_t q = 0; q < queries; ++q) {
            extract_query(positions[q]);
            checksum += output[0];
        }
        const auto end_time = BenchmarkClock::now();
        times.push_back(1000.0 * elapsed_ms(begin_time, end_time) / queries);
        for (size_t q = 0; q < std::min<size_t>(queries, 100); ++q) {
            extract_query(positions[q]);
            if (!std::equal(output.begin(), output.end(),
                            original + positions[q]))
                throw std::runtime_error("Block extraction correctness check failed");
        }
        benchmark_sink = checksum;
    }
    return times;
}

template<class RankStructure, class Extractor>
std::vector<double> benchmark_rlz_series(RankStructure& rank,
                                         const uint8_t* text, size_t n,
                                         Extractor extract) {
    constexpr int queries = 10000;
    std::vector<double> times;
    std::mt19937 rng(12345);
    std::vector<int> positions(queries);
    std::vector<uint8_t> buffer;
    for (int length = 8; length <= 4096; length *= 2) {
        if (static_cast<size_t>(length) > n) { times.push_back(0.0); continue; }
        buffer.resize(length);
        std::uniform_int_distribution<int> distribution(0, static_cast<int>(n) - length);
        for (int& position : positions) position = distribution(rng);
        uint64_t checksum = 0;
        const auto begin = BenchmarkClock::now();
        for (int position : positions) {
            const int copied = extract(rank, position, position + length - 1,
                                       buffer.data());
            checksum += copied > 0 ? buffer[0] : 0;
        }
        const auto end = BenchmarkClock::now();
        times.push_back(1000.0 * elapsed_ms(begin, end) / queries);
        benchmark_sink = checksum;
    }
    return times;
}

int run_figure1_benchmark(const std::string& list_name) {
    const auto inputs = read_size_benchmark_inputs(list_name);
    const std::filesystem::path output_name =
        std::filesystem::absolute(list_name).parent_path() / "figure1_table.csv";
    std::ofstream csv(output_name);
    if (!csv) throw std::runtime_error("Cannot create " + output_name.string());
    csv << "File,Method,Block/frame size,Compressed size (bytes),Compression ratio";
    for (int length = 8; length <= 4096; length *= 2) csv << ',' << length;
    csv << '\n' << std::setprecision(10);
    ZstdApi zstd;

    const auto write_row = [&](const std::string& file, const std::string& method,
                               size_t parameter, uint64_t compressed, size_t original,
                               const std::vector<double>& times) {
        csv << csv_field(std::filesystem::path(file).filename().string()) << ','
            << csv_field(method) << ',' << parameter << ',' << compressed << ','
            << static_cast<double>(original) / compressed;
        for (double value : times) csv << ',' << value;
        csv << '\n';
    };

    for (const SizeBenchmarkInput& input : inputs) {
        std::cout << "Figure 1 benchmark: " << input.file << std::endl;
        WordBasedText text(input.file, 1);
        const size_t n = static_cast<size_t>(text.Nchar);
        size_t dictionary_size = 1;
        while (static_cast<long double>(dictionary_size) <= 0.014L * n)
            dictionary_size *= 2;
        std::cout << "  text=" << n << " dictionary=" << dictionary_size
                  << " L=" << input.max_suffix_length << std::endl;
        TextParser parser(&text);
        auto stage_begin = BenchmarkClock::now();
        std::cout << "  Building FM-index dictionary..." << std::flush;
        parser.build_suffix_array_dictionary(dictionary_size,
                                             input.max_suffix_length, -1,
                                             input.seed);
        std::cout << " done ("
                  << elapsed_ms(stage_begin, BenchmarkClock::now())
                  << " ms)" << std::endl;
        std::cout << "  Constructing RLZ parse..." << std::flush;
        stage_begin = BenchmarkClock::now();
        parser.process_text_liao(1, false);
        std::cout << " done ("
                  << elapsed_ms(stage_begin, BenchmarkClock::now())
                  << " ms, " << parser.segments_number() << " phrases)"
                  << std::endl;
        RLZ_Rank bit_fixed(parser.lengths_data(), parser.offsets_data(),
                           parser.dictionary_data(), parser.segments_number(),
                           text.Nchar, parser.dictionary_size(), false,
                           0, 0, 1, false);
        RLZ_Rank_ef ef_vbyte(parser.lengths_data(), parser.offsets_data(),
                            parser.dictionary_data(), parser.segments_number(),
                            text.Nchar, parser.dictionary_size(), 0, 0, 1, true);
        const uint64_t bit_size = dictionary_size + bit_fixed.space_usage() +
                                  bit_fixed.packed_offsets_size();
        const uint64_t ef_size = dictionary_size + ef_vbyte.space_usage() +
                                 ef_vbyte.packed_offsets_size();
        std::cout << "  Timing Bit-Fix queries (10,000 per length)..."
                  << std::flush;
        const auto bit_times = benchmark_rlz_series(
            bit_fixed, parser.text_data(), n,
            [](RLZ_Rank& r, int s, int e, uint8_t* out) {
                return r.get_seg_word_scan(s, e, out);
            });
        std::cout << " done" << std::endl;
        std::cout << "  Timing EF-VB queries (10,000 per length)..."
                  << std::flush;
        const auto ef_times = benchmark_rlz_series(
            ef_vbyte, parser.text_data(), n,
            [](RLZ_Rank_ef& r, int s, int e, uint8_t* out) {
                return r.get_seg(s, e, out);
            });
        std::cout << " done" << std::endl;
        write_row(input.file, "Our Bit-Fix", 0, bit_size, n, bit_times);
        write_row(input.file, "Our EF-VB", 0, ef_size, n, ef_times);

        for (const auto& target : std::array<std::pair<const char*, uint64_t>, 2>{{
                 {"Bit-Fix", bit_size}, {"EF-VB", ef_size}}}) {
            std::cout << "  Matching Zstd Seekable to " << target.first
                      << " size..." << std::flush;
            ZstdBlockArchive seekable = closest_seekable_archive(
                zstd, parser.text_data(), n, target.second);
            std::cout << " frame=" << seekable.block_size
                      << " bytes, error="
                      << 100.0 * std::abs(
                             static_cast<double>(seekable.size()) -
                             static_cast<double>(target.second)) /
                             static_cast<double>(target.second)
                      << "%; timing queries..." << std::flush;
            write_row(input.file, std::string("Zstd Seekable ") + target.first,
                      seekable.block_size, seekable.size(), n,
                      benchmark_block_archive(zstd, seekable, parser.text_data()));
            std::cout << " done" << std::endl;
            std::cout << "  Selecting SquashFS block for " << target.first
                      << "..." << std::flush;
            ZstdBlockArchive squash = closest_squashfs_archive(
                zstd, parser.text_data(), n, target.second,
                seekable.block_size);
            std::cout << " block=" << squash.block_size
                      << " bytes; timing queries..." << std::flush;
            write_row(input.file, std::string("SquashFS ") + target.first,
                      squash.block_size, squash.size(), n,
                      benchmark_block_archive(zstd, squash, parser.text_data()));
            std::cout << " done" << std::endl;
        }
        csv.flush();
    }
    std::cout << "Figure 1 table: " << output_name.string() << std::endl;
    return 0;
}

} // namespace

int main(int argc, char** argv) {
    if ((argc != 2 && argc != 3) ||
        std::filesystem::path(argv[1]).extension() != ".csv") {
        std::cerr << "Usage: rlz_benchmark <files.csv> "
                     "[table1|table2|figure1]\n";
        return 1;
    }

    try {
        if (argc == 2) {
            run_dictionary_quality_benchmark(argv[1]);
            run_structure_size_benchmark(argv[1]);
            return run_figure1_benchmark(argv[1]);
        }

        const std::string mode = argv[2];
        if (mode == "table1")
            return run_dictionary_quality_benchmark(argv[1]);
        if (mode == "table2")
            return run_structure_size_benchmark(argv[1]);
        if (mode == "figure1")
            return run_figure1_benchmark(argv[1]);
        throw std::invalid_argument(
            "Unknown mode '" + mode + "'. Expected table1, table2, or figure1.");
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
