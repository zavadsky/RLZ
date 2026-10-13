#include <iostream>
#include <algorithm>
#include <cstring>
#include <stdexcept>
#define PASTA_SELECT_IMPLEMENTED_IN_HEADER
#include "include/Rank.hpp"

using namespace std;

namespace {

int previous_one_wordwise(const pasta::BitVector& bv, int position) {
    size_t word_index = static_cast<size_t>(position) >> 6;
    const unsigned bit_index = static_cast<unsigned>(position) & 63U;
    uint64_t word = bv.data(word_index);
    word &= bit_index == 63U ? ~uint64_t{0}
                             : ((uint64_t{1} << (bit_index + 1U)) - 1U);

    while (word == 0) {
        if (word_index == 0)
            return -1;
        word = bv.data(--word_index);
    }
    return static_cast<int>((word_index << 6) +
                            (63U - static_cast<unsigned>(__builtin_clzll(word))));
}

int next_one_wordwise(const pasta::BitVector& bv, int position,
                      int text_length) {
    if (position >= text_length)
        return text_length;

    size_t word_index = static_cast<size_t>(position) >> 6;
    const unsigned bit_index = static_cast<unsigned>(position) & 63U;
    uint64_t word = bv.data(word_index) & (~uint64_t{0} << bit_index);
    const size_t word_count = (static_cast<size_t>(text_length) + 63U) >> 6;

    while (word == 0) {
        if (++word_index >= word_count)
            return text_length;
        word = bv.data(word_index);
    }

    const size_t result = (word_index << 6) +
                          static_cast<unsigned>(__builtin_ctzll(word));
    return result < static_cast<size_t>(text_length)
               ? static_cast<int>(result)
               : text_length;
}

inline void copy_with_64bit_chunks(uint8_t* destination,
                                   const uint8_t* source,
                                   size_t count) {
    if (count < sizeof(uint64_t)) {
        std::memcpy(destination, source, count);
        return;
    }

    size_t position = 0;
    while (position + sizeof(uint64_t) <= count) {
        uint64_t value;
        std::memcpy(&value, source + position, sizeof(value));
        std::memcpy(destination + position, &value, sizeof(value));
        position += sizeof(uint64_t);
    }

    // Cover a 1--7 byte tail with one overlapping 8-byte copy.  Both writes
    // use the same source positions for the overlap, so the result is exact.
    if (position < count) {
        uint64_t value;
        const size_t tail_begin = count - sizeof(uint64_t);
        std::memcpy(&value, source + tail_begin, sizeof(value));
        std::memcpy(destination + tail_begin, &value, sizeof(value));
    }
}

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

} // namespace

RLZ_Rank::RLZ_Rank(const uint32_t* lengths, const uint32_t* offsets,
                   const uint8_t* dictionary, int segments, int text_length,
                   uint64_t dictionary_size, bool verbose,
                   uint64_t common_size,
                   uint64_t regular_block_part_size,
                   int dictionary_blocks, bool use_vbyte_offsets)
    : dictionary(dictionary),
      segments(segments), text_length(text_length), common_size(common_size),
      regular_block_part_size(regular_block_part_size),
      dictionary_blocks(dictionary_blocks), vbyte_offsets(use_vbyte_offsets)
{
    uint64_t logical_dictionary_size = dictionary_size;
    if (common_size > 0 && dictionary_blocks > 1) {
        const uint64_t last_part = dictionary_size - common_size -
            regular_block_part_size * static_cast<uint64_t>(dictionary_blocks - 1);
        logical_dictionary_size = common_size +
            std::max(regular_block_part_size, last_part);
    }
    if (vbyte_offsets)
        packed_offsets.build(offsets, segments, logical_dictionary_size);
    else
        fixed_offsets.build(offsets, segments, logical_dictionary_size);

    bv.resize(text_length, 0);
    if (verbose)
        cout << "N=" << text_length << "; phrases=" << segments << endl;

    // A 1 marks the beginning of every phrase.
    size_t phrase_begin = 0;
    for (int phrase = 0; phrase < segments; ++phrase) {
        if (phrase_begin >= static_cast<size_t>(text_length))
            break;
        bv[phrase_begin] = 1;
        phrase_begin += lengths[phrase];
    }
    bvr = new pasta::Rank<pasta::OptimizedFor::ONE_QUERIES>(bv);
    if (verbose) {
        cout << "Boundary bitvector size=" << bv.space_usage() << " bytes\n"
             << "Rank support size=" << bvr->space_usage() << " bytes\n";
    }
}

RLZ_Rank::~RLZ_Rank()
{
    delete bvr;
}

uint64_t RLZ_Rank::space_usage() const {
    return static_cast<uint64_t>(bv.space_usage())
         + static_cast<uint64_t>(bvr->space_usage());
}

uint64_t RLZ_Rank::packed_offsets_size() const {
    return vbyte_offsets ? packed_offsets.space_usage()
                         : fixed_offsets.space_usage();
}

uint64_t RLZ_Rank::offset_low_size() const {
    return vbyte_offsets ? packed_offsets.low_size() : fixed_offsets.space_usage();
}

uint64_t RLZ_Rank::offset_presence_size() const {
    return vbyte_offsets ? packed_offsets.presence_size() : 0;
}

uint64_t RLZ_Rank::offset_ternary_size() const {
    return vbyte_offsets ? packed_offsets.ternary_size() : 0;
}

unsigned RLZ_Rank::offset_width() const {
    return vbyte_offsets ? packed_offsets.width() : fixed_offsets.width();
}

uint64_t RLZ_Rank::offset_at(int phrase) const {
    return vbyte_offsets ? packed_offsets.at(phrase) : fixed_offsets.at(phrase);
}

uint64_t RLZ_Rank::physical_offset(uint64_t logical_offset,
                                   int phrase_begin) const {
    if (common_size == 0 || dictionary_blocks <= 1 ||
        logical_offset < common_size)
        return logical_offset;
    const uint64_t text_block_size =
        (static_cast<uint64_t>(text_length) + dictionary_blocks - 1) /
        static_cast<uint64_t>(dictionary_blocks);
    const uint64_t block = std::min<uint64_t>(
        static_cast<uint64_t>(dictionary_blocks - 1),
        static_cast<uint64_t>(phrase_begin) / text_block_size);
    return logical_offset + block * regular_block_part_size;
}

int RLZ_Rank::get_seg(int start, int end, uint8_t* buffer) const {
    if (buffer == nullptr || text_length <= 0 || start < 0 || start >= text_length || end < start)
        return 0;

    end = std::min(end, text_length - 1);

    // rank is used only once, to identify the first phrase.
    int phrase = static_cast<int>(bvr->rank1(start + 1)) - 1;

    int phrase_begin = start;
    while (phrase_begin > 0 && !bv[phrase_begin])
        --phrase_begin;

    int phrase_end = phrase_begin + 1;
    while (phrase_end < text_length && !bv[phrase_end])
        ++phrase_end;

    int copied = 0;
    int current = start;

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

        // Subsequent phrases are found only by scanning to the next 1-bit.
        ++phrase;
        phrase_begin = phrase_end;
        phrase_end = phrase_begin + 1;
        while (phrase_end < text_length && !bv[phrase_end])
            ++phrase_end;
    }

    return copied;
}

int RLZ_Rank::get_seg_word_scan(int start, int end, uint8_t* buffer) const {
    if (buffer == nullptr || text_length <= 0 || start < 0 ||
        start >= text_length || end < start)
        return 0;

    end = std::min(end, text_length - 1);
    int phrase = static_cast<int>(bvr->rank1(start + 1)) - 1;
    if (phrase < 0 || phrase >= segments)
        return 0;

    int phrase_begin = previous_one_wordwise(bv, start);
    int phrase_end = next_one_wordwise(bv, start + 1, text_length);
    if (phrase_begin < 0 || phrase_end <= start)
        return 0;

    int copied = 0;
    int current = start;
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
        phrase_end = next_one_wordwise(bv, phrase_begin + 1, text_length);
        if (phrase_end <= phrase_begin)
            return 0;
    }
    return copied;
}

int RLZ_Rank::get_seg_word_scan_copy8(int start, int end,
                                      uint8_t* buffer) const {
    if (buffer == nullptr || text_length <= 0 || start < 0 ||
        start >= text_length || end < start)
        return 0;

    end = std::min(end, text_length - 1);
    int phrase = static_cast<int>(bvr->rank1(start + 1)) - 1;
    if (phrase < 0 || phrase >= segments)
        return 0;

    int phrase_begin = previous_one_wordwise(bv, start);
    int phrase_end = next_one_wordwise(bv, start + 1, text_length);
    if (phrase_begin < 0 || phrase_end <= start)
        return 0;

    int copied = 0;
    int current = start;
    while (current <= end && phrase < segments) {
        const int within_phrase = current - phrase_begin;
        const int copy_end = std::min(end + 1, phrase_end);
        const int count = copy_end - current;

        if (common_size > 0 && dictionary_blocks > 1) {
            copy_logical_dictionary(buffer + copied, dictionary,
                offset_at(phrase) + within_phrase, static_cast<size_t>(count),
                phrase_begin, text_length, common_size,
                regular_block_part_size, dictionary_blocks);
        } else {
            copy_with_64bit_chunks(buffer + copied,
                                   dictionary + offset_at(phrase) + within_phrase,
                                   static_cast<size_t>(count));
        }
        copied += count;
        current = copy_end;
        if (current > end)
            break;

        ++phrase;
        phrase_begin = phrase_end;
        phrase_end = next_one_wordwise(bv, phrase_begin + 1, text_length);
        if (phrase_end <= phrase_begin)
            return 0;
    }
    return copied;
}

bool RLZ_Rank::check_segment(int start, int end, const uint8_t* true_text,
                             uint8_t* buffer) const {
    if (true_text == nullptr || buffer == nullptr || start < 0 ||
        start >= text_length || end < start) {
        std::cerr << "Invalid segment range [" << start << ", " << end << "]"
                  << std::endl;
        return false;
    }

    end = std::min(end, text_length - 1);
    const int expected_length = end - start + 1;
    const int copied = get_seg(start, end, buffer);

    if (copied != expected_length) {
        std::cerr << "Incorrect extracted length for [" << start << ", " << end
                  << "]: expected " << expected_length << ", got " << copied
                  << std::endl;
        return false;
    }

    for (int i = 0; i < expected_length; ++i) {
        if (buffer[i] != true_text[start + i]) {
            std::cerr << "Mismatch for [" << start << ", " << end
                      << "] at text position " << start + i
                      << ": expected " << static_cast<unsigned>(true_text[start + i])
                      << ", got " << static_cast<unsigned>(buffer[i]) << std::endl;
            return false;
        }
    }

    return true;
}

bool RLZ_Rank::check_segment_word_scan(int start, int end,
                                       const uint8_t* true_text,
                                       uint8_t* buffer) const {
    if (true_text == nullptr || buffer == nullptr || start < 0 ||
        start >= text_length || end < start)
        return false;

    end = std::min(end, text_length - 1);
    const int expected_length = end - start + 1;
    const int copied = get_seg_word_scan(start, end, buffer);
    if (copied != expected_length)
        return false;

    for (int i = 0; i < expected_length; ++i) {
        if (buffer[i] != true_text[start + i]) {
            std::cerr << "Word-scan mismatch for [" << start << ", " << end
                      << "] at text position " << start + i << std::endl;
            return false;
        }
    }
    return true;
}

bool RLZ_Rank::check_segment_word_scan_copy8(int start, int end,
                                             const uint8_t* true_text,
                                             uint8_t* buffer) const {
    if (true_text == nullptr || buffer == nullptr || start < 0 ||
        start >= text_length || end < start)
        return false;

    end = std::min(end, text_length - 1);
    const int expected_length = end - start + 1;
    const int copied = get_seg_word_scan_copy8(start, end, buffer);
    if (copied != expected_length)
        return false;

    for (int i = 0; i < expected_length; ++i) {
        if (buffer[i] != true_text[start + i]) {
            std::cerr << "64-bit-copy mismatch for [" << start << ", " << end
                      << "] at text position " << start + i << std::endl;
            return false;
        }
    }
    return true;
}

#ifndef PASTA_SELECT_IMPLEMENTED_IN_HEADER
namespace pasta {
constexpr uint8_t kSelectInByte[2048] = {
    8, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3,
    0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0,
    1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1,
    0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0,
    2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2,
    0, 1, 0, 7, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0,
    1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1,
    0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 6, 0, 1, 0, 2, 0, 1, 0,
    3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 5,
    0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0, 4, 0, 1, 0, 2, 0, 1, 0, 3, 0,
    1, 0, 2, 0, 1, 0, 8, 8, 8, 1, 8, 2, 2, 1, 8, 3, 3, 1, 3, 2, 2, 1, 8, 4, 4,
    1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1,
    3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 6, 6, 1, 6,
    2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2,
    2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2,
    1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 7, 7, 1, 7, 2, 2, 1, 7, 3, 3, 1, 3, 2, 2, 1,
    7, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 5, 5, 1, 5, 2, 2, 1, 5,
    3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 7, 6,
    6, 1, 6, 2, 2, 1, 6, 3, 3, 1, 3, 2, 2, 1, 6, 4, 4, 1, 4, 2, 2, 1, 4, 3, 3,
    1, 3, 2, 2, 1, 6, 5, 5, 1, 5, 2, 2, 1, 5, 3, 3, 1, 3, 2, 2, 1, 5, 4, 4, 1,
    4, 2, 2, 1, 4, 3, 3, 1, 3, 2, 2, 1, 8, 8, 8, 8, 8, 8, 8, 2, 8, 8, 8, 3, 8,
    3, 3, 2, 8, 8, 8, 4, 8, 4, 4, 2, 8, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 5, 8, 5,
    5, 2, 8, 5, 5, 3, 5, 3, 3, 2, 8, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3,
    2, 8, 8, 8, 6, 8, 6, 6, 2, 8, 6, 6, 3, 6, 3, 3, 2, 8, 6, 6, 4, 6, 4, 4, 2,
    6, 4, 4, 3, 4, 3, 3, 2, 8, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3, 3, 2, 6,
    5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 7, 8, 7, 7, 2, 8, 7,
    7, 3, 7, 3, 3, 2, 8, 7, 7, 4, 7, 4, 4, 2, 7, 4, 4, 3, 4, 3, 3, 2, 8, 7, 7,
    5, 7, 5, 5, 2, 7, 5, 5, 3, 5, 3, 3, 2, 7, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3,
    4, 3, 3, 2, 8, 7, 7, 6, 7, 6, 6, 2, 7, 6, 6, 3, 6, 3, 3, 2, 7, 6, 6, 4, 6,
    4, 4, 2, 6, 4, 4, 3, 4, 3, 3, 2, 7, 6, 6, 5, 6, 5, 5, 2, 6, 5, 5, 3, 5, 3,
    3, 2, 6, 5, 5, 4, 5, 4, 4, 2, 5, 4, 4, 3, 4, 3, 3, 2, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 3, 8, 8, 8, 8, 8, 8, 8, 4, 8, 8, 8, 4, 8, 4, 4, 3,
    8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 3, 8, 8, 8, 5, 8, 5, 5, 4, 8,
    5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 3, 8, 8,
    8, 6, 8, 6, 6, 4, 8, 6, 6, 4, 6, 4, 4, 3, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6,
    5, 6, 5, 5, 3, 8, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 8,
    8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 3, 8, 8, 8, 7, 8, 7, 7, 4, 8, 7, 7, 4, 7,
    4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 3, 8, 7, 7, 5, 7, 5,
    5, 4, 7, 5, 5, 4, 5, 4, 4, 3, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6,
    3, 8, 7, 7, 6, 7, 6, 6, 4, 7, 6, 6, 4, 6, 4, 4, 3, 8, 7, 7, 6, 7, 6, 6, 5,
    7, 6, 6, 5, 6, 5, 5, 3, 7, 6, 6, 5, 6, 5, 5, 4, 6, 5, 5, 4, 5, 4, 4, 3, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8,
    8, 8, 8, 8, 5, 8, 8, 8, 5, 8, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 4, 8, 8, 8, 8, 8,
    8, 8, 6, 8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 6, 8, 6, 6, 5, 8, 6, 6, 5, 6, 5,
    5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8,
    7, 8, 8, 8, 7, 8, 7, 7, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5,
    8, 8, 8, 7, 8, 7, 7, 5, 8, 7, 7, 5, 7, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 7, 8,
    8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 4, 8, 8,
    8, 7, 8, 7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 7, 7, 6, 7, 6, 6, 5, 7, 6, 6,
    5, 6, 5, 5, 4, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 5, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 6,
    8, 8, 8, 6, 8, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 5, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7,
    8, 7, 7, 6, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 7, 8,
    7, 7, 6, 8, 7, 7, 6, 7, 6, 6, 5, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 6, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    7, 8, 8, 8, 8, 8, 8, 8, 7, 8, 8, 8, 7, 8, 7, 7, 6, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
    8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 7};

/*[[nodiscard]]*/ uint64_t select_pasta1(uint64_t x, uint64_t k) {
#if !defined(__BMI2__)
    // Portable zero-based select within one 64-bit word.
    if (x == 0 || k >= static_cast<uint64_t>(__builtin_popcountll(x)))
        return 64;
    while (k-- != 0)
        x &= x - 1;
    return static_cast<uint64_t>(__builtin_ctzll(x));
#elif defined(__GNUC__) || defined(__clang__)

  // GCC and Clang won't inline the intrinsics.
  uint64_t result = uint64_t(1) << k;

  asm("pdep %1, %0, %0\n\t"
      "tzcnt %0, %0"
      : "+r"(result)
      : "r"(x));

  return result;
#else
  return _tzcnt_u64(_pdep_u64(1ULL << k, x));
#endif
}
}
#endif
