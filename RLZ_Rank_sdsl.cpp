#include "include/RLZ_Rank_sdsl.hpp"
#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/rank_select.hpp"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <stdexcept>

// Compatibility helper for PASTA releases in which RankSelect declares
// select_pasta1() in the header but leaves its definition to the application.
// k is a zero-based rank within the set bits of x.
namespace pasta {
uint64_t select_pasta1(uint64_t x, uint64_t k) {
    if (x == 0 || k >= static_cast<uint64_t>(__builtin_popcountll(x)))
        return 64;
    while (k-- != 0)
        x &= x - 1;
    return static_cast<uint64_t>(__builtin_ctzll(x));
}
} // namespace pasta

using namespace std;

struct RLZ_Rank_ef::EFData {
    using HighSelect =
        pasta::RankSelect<pasta::OptimizedFor::ZERO_QUERIES>;

    pasta::BitVector high;
    HighSelect* high_select = nullptr;
    std::vector<uint64_t> low;
    unsigned low_bits = 0;

    ~EFData() {
        delete high_select;
    }
};

namespace {

void copy_logical_dictionary(uint8_t* destination, const uint8_t* dictionary,
                             uint64_t logical_offset, size_t count,
                             int phrase_begin, int text_length,
                             uint64_t common_size,
                             uint64_t regular_part_size, int blocks) {
    if (common_size == 0 || blocks <= 1) {
        std::memcpy(destination, dictionary + logical_offset, count);
        return;
    }
    const uint64_t block_length =
        (static_cast<uint64_t>(text_length) + blocks - 1) /
        static_cast<uint64_t>(blocks);
    const uint64_t block = std::min<uint64_t>(blocks - 1,
        static_cast<uint64_t>(phrase_begin) / block_length);
    size_t copied = 0;
    if (logical_offset < common_size) {
        const size_t prefix = std::min<size_t>(
            count, static_cast<size_t>(common_size - logical_offset));
        std::memcpy(destination, dictionary + logical_offset, prefix);
        logical_offset += prefix;
        copied += prefix;
    }
    if (copied < count) {
        const uint64_t physical = logical_offset + block * regular_part_size;
        std::memcpy(destination + copied, dictionary + physical,
                    count - copied);
    }
}

uint64_t previous_zero_wordwise(const pasta::BitVector& bv,
                                uint64_t position) {
    uint64_t word_index = position >> 6;
    const unsigned bit_index = static_cast<unsigned>(position & 63);
    uint64_t word = ~bv.data(static_cast<size_t>(word_index));
    word &= bit_index == 63U ? ~uint64_t{0}
                             : ((uint64_t{1} << (bit_index + 1U)) - 1U);

    while (word == 0) {
        if (word_index == 0)
            return uint64_t(-1);
        word = ~bv.data(static_cast<size_t>(--word_index));
    }
    return (word_index << 6) +
           (63U - static_cast<unsigned>(__builtin_clzll(word)));
}

} // namespace


RLZ_Rank_ef::RLZ_Rank_ef(const uint32_t* lengths,
                         const uint32_t* offsets,
                         const uint8_t* dictionary, int segments,
                         int text_length, uint64_t dictionary_size,
                         uint64_t common_size,
                         uint64_t regular_block_part_size,
                         int dictionary_blocks, bool use_vbyte_offsets)
    : ef_data(new EFData), dictionary(dictionary),
      segments(segments), text_length(text_length), common_size(common_size),
      regular_block_part_size(regular_block_part_size),
      dictionary_blocks(dictionary_blocks), vbyte_offsets(use_vbyte_offsets)
{
    uint64_t logical_size = dictionary_size;
    if (common_size > 0 && dictionary_blocks > 1) {
        const uint64_t last = dictionary_size - common_size -
            regular_block_part_size * static_cast<uint64_t>(dictionary_blocks - 1);
        logical_size = common_size + std::max(regular_block_part_size, last);
    }
    if (vbyte_offsets)
        packed_offsets.build(offsets, segments, logical_size);
    else
        fixed_offsets.build(offsets, segments, logical_size);

    const uint64_t average_gap = segments > 0
        ? static_cast<uint64_t>(text_length) /
              static_cast<uint64_t>(segments)
        : 1;
    while (ef_data->low_bits < 63 &&
           (uint64_t{1} << (ef_data->low_bits + 1)) <= average_gap)
        ++ef_data->low_bits;

    const uint64_t low_bit_count =
        static_cast<uint64_t>(segments) * ef_data->low_bits;
    ef_data->low.assign((low_bit_count + 63) / 64, 0);
    const uint64_t low_mask = ef_data->low_bits == 0
        ? 0
        : ((uint64_t{1} << ef_data->low_bits) - 1);

    const uint64_t high_universe =
        (static_cast<uint64_t>(text_length) >> ef_data->low_bits) +
        static_cast<uint64_t>(segments) + 1;
    ef_data->high.resize(static_cast<size_t>(high_universe), 0);

    uint64_t phrase_begin = 0;
    for (int phrase = 0; phrase < segments; ++phrase) {
        const uint64_t high = phrase_begin >> ef_data->low_bits;
        ef_data->high[static_cast<size_t>(high + phrase)] = 1;

        if (ef_data->low_bits != 0) {
            const uint64_t low = phrase_begin & low_mask;
            const uint64_t bit_position =
                static_cast<uint64_t>(phrase) * ef_data->low_bits;
            const size_t word = static_cast<size_t>(bit_position >> 6);
            const unsigned shift = static_cast<unsigned>(bit_position & 63);
            ef_data->low[word] |= low << shift;
            if (shift + ef_data->low_bits > 64)
                ef_data->low[word + 1] |= low >> (64 - shift);
        }
        phrase_begin += lengths[phrase];
    }
    ef_data->high_select = new EFData::HighSelect(ef_data->high);
}

RLZ_Rank_ef::~RLZ_Rank_ef() {
    delete ef_data;
}

uint64_t RLZ_Rank_ef::offset_at(int phrase) const {
    return vbyte_offsets ? packed_offsets.at(phrase) : fixed_offsets.at(phrase);
}

uint64_t RLZ_Rank_ef::physical_offset(uint64_t logical_offset,
                                      int phrase_begin) const {
    if (common_size == 0 || dictionary_blocks <= 1 || logical_offset < common_size)
        return logical_offset;
    const uint64_t block_length =
        (static_cast<uint64_t>(text_length) + dictionary_blocks - 1) /
        static_cast<uint64_t>(dictionary_blocks);
    const uint64_t block = std::min<uint64_t>(dictionary_blocks - 1,
        static_cast<uint64_t>(phrase_begin) / block_length);
    return logical_offset + block * regular_block_part_size;
}

uint64_t RLZ_Rank_ef::low_at(int phrase) const {
    if (ef_data->low_bits == 0)
        return 0;
    const uint64_t bit_position =
        static_cast<uint64_t>(phrase) * ef_data->low_bits;
    const size_t word = static_cast<size_t>(bit_position >> 6);
    const unsigned shift = static_cast<unsigned>(bit_position & 63);
    uint64_t value = ef_data->low[word] >> shift;
    if (shift + ef_data->low_bits > 64)
        value |= ef_data->low[word + 1] << (64 - shift);
    return value & ((uint64_t{1} << ef_data->low_bits) - 1);
}

uint64_t RLZ_Rank_ef::value_at(int phrase, uint64_t high_position) const {
    const uint64_t high = high_position - static_cast<uint64_t>(phrase);
    return (high << ef_data->low_bits) | low_at(phrase);
}

uint64_t RLZ_Rank_ef::next_high_one(uint64_t position) const {
    const uint64_t high_size = ef_data->high.size();
    if (position >= high_size)
        return high_size;

    uint64_t word_index = position >> 6;
    const unsigned bit_index = static_cast<unsigned>(position & 63);
    uint64_t word = ef_data->high.data(static_cast<size_t>(word_index)) &
                    (~uint64_t{0} << bit_index);
    const uint64_t word_count = (high_size + 63) >> 6;

    while (word == 0) {
        if (++word_index >= word_count)
            return high_size;
        word = ef_data->high.data(static_cast<size_t>(word_index));
    }
    const uint64_t result = (word_index << 6) +
                            static_cast<unsigned>(__builtin_ctzll(word));
    return std::min(result, high_size);
}

int RLZ_Rank_ef::predecessor(int position,
                             uint64_t& high_position) const {
    // Count all phrase beginnings strictly smaller than position + 1.
    const uint64_t x = static_cast<uint64_t>(position) + 1;
    const uint64_t high = x >> ef_data->low_bits;
    const uint64_t low_mask = ef_data->low_bits == 0
        ? 0
        : ((uint64_t{1} << ef_data->low_bits) - 1);
    const uint64_t low = x & low_mask;

    // The (high+1)-st zero terminates the bucket for this high part.
    const uint64_t bucket_end_position =
        ef_data->high_select->select0(high + 1);
    const uint64_t bucket_end = bucket_end_position - high;

    // Locate the preceding zero locally; the intervening one-run is exactly
    // the current Elias--Fano bucket.
    const uint64_t previous_zero = high == 0
        ? uint64_t(-1)
        : previous_zero_wordwise(ef_data->high,
                                 bucket_end_position - 1);
    const uint64_t bucket_begin = high == 0
        ? 0
        : previous_zero - (high - 1);

    // Find the first value in the bucket whose low part is >= low.
    uint64_t left = bucket_begin;
    uint64_t right = bucket_end;
    while (left < right) {
        const uint64_t middle = left + ((right - left) >> 1);
        if (low_at(static_cast<int>(middle)) < low)
            left = middle + 1;
        else
            right = middle;
    }

    const int phrase = static_cast<int>(left) - 1;
    if (phrase < 0)
        return -1;

    if (static_cast<uint64_t>(phrase) >= bucket_begin) {
        high_position = high + static_cast<uint64_t>(phrase);
    } else {
        high_position = next_high_one(0);
        if (previous_zero != uint64_t(-1) && previous_zero > 0) {
            // The predecessor belongs to an earlier bucket.  Scan left from
            // the zero preceding the current bucket to its last one-bit.
            uint64_t word_index = (previous_zero - 1) >> 6;
            const unsigned bit_index =
                static_cast<unsigned>((previous_zero - 1) & 63);
            uint64_t word =
                ef_data->high.data(static_cast<size_t>(word_index));
            word &= bit_index == 63U
                ? ~uint64_t{0}
                : ((uint64_t{1} << (bit_index + 1U)) - 1U);
            while (word == 0) {
                word = ef_data->high.data(static_cast<size_t>(--word_index));
            }
            high_position = (word_index << 6) +
                (63U - static_cast<unsigned>(__builtin_clzll(word)));
        }
    }
    return phrase;
}

int RLZ_Rank_ef::get_seg(int start, int end, uint8_t* buffer) const {
    if (buffer == nullptr || text_length <= 0 || start < 0 ||
        start >= text_length || end < start)
        return 0;

    end = std::min(end, text_length - 1);

    // One zero-select followed by a binary search initializes the EF cursor.
    uint64_t high_position = 0;
    int phrase = predecessor(start, high_position);
    if (phrase < 0 || phrase >= segments)
        return 0;
    int phrase_begin = static_cast<int>(value_at(phrase, high_position));

    uint64_t next_high_position = ef_data->high.size();
    int phrase_end = text_length;
    if (phrase + 1 < segments) {
        next_high_position = next_high_one(high_position + 1);
        phrase_end = static_cast<int>(value_at(phrase + 1,
                                               next_high_position));
    }

    int current = start;
    int copied = 0;
    while (current <= end && phrase < segments) {
        const int within_phrase = current - phrase_begin;
        const int copy_end = std::min(end + 1, phrase_end);
        const int count = copy_end - current;
        copy_logical_dictionary(buffer + copied, dictionary,
            offset_at(phrase) + within_phrase, static_cast<size_t>(count),
            phrase_begin, text_length, common_size, regular_block_part_size,
            dictionary_blocks);
        copied += count;
        current = copy_end;

        if (current > end)
            break;

        ++phrase;
        phrase_begin = phrase_end;
        high_position = next_high_position;
        if (phrase + 1 < segments) {
            next_high_position = next_high_one(high_position + 1);
            phrase_end = static_cast<int>(value_at(phrase + 1,
                                                   next_high_position));
        } else {
            phrase_end = text_length;
        }
    }
    return copied;
}

bool RLZ_Rank_ef::check_segment(int start, int end,
                                const uint8_t* true_text,
                                uint8_t* buffer) const {
    if (true_text == nullptr || buffer == nullptr || start < 0 ||
        start >= text_length || end < start)
        return false;
    end = std::min(end, text_length - 1);
    const int expected = end - start + 1;
    const int copied = get_seg(start, end, buffer);
    return copied == expected &&
           std::memcmp(buffer, true_text + start,
                       static_cast<size_t>(expected)) == 0;
}

uint64_t RLZ_Rank_ef::space_usage() const {
    using RankBase =
        pasta::Rank<pasta::OptimizedFor::ZERO_QUERIES, pasta::BitVector>;
    const uint64_t high_support_size =
        static_cast<uint64_t>(ef_data->high_select->space_usage()) +
        static_cast<uint64_t>(
            ef_data->high_select->RankBase::space_usage() - sizeof(RankBase));
    return static_cast<uint64_t>(ef_data->high.space_usage()) +
           static_cast<uint64_t>(ef_data->low.size()) * sizeof(uint64_t) +
           high_support_size;
}

uint64_t RLZ_Rank_ef::packed_offsets_size() const {
    return vbyte_offsets ? packed_offsets.space_usage()
                         : fixed_offsets.space_usage();
}

unsigned RLZ_Rank_ef::offset_width() const {
    return vbyte_offsets ? packed_offsets.width() : fixed_offsets.width();
}



