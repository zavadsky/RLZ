#pragma once

#ifndef RANK_HPP
#define RANK_HPP
#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/rank.hpp"
#include "PackedTernaryOffsets.hpp"
#include "PackedFixedOffsets.hpp"
#include <cstdint>
#include <vector>


class RLZ_Rank
{
    public:
        RLZ_Rank(const uint32_t* lengths, const uint32_t* offsets,
                 const uint8_t* dictionary, int segments, int text_length,
                 uint64_t dictionary_size, bool verbose = true,
                 uint64_t common_size = 0,
                 uint64_t regular_block_part_size = 0,
                 int dictionary_blocks = 1,
                 bool vbyte_offsets = true);
        virtual ~RLZ_Rank();
        int get_seg(int start, int end, uint8_t* buffer) const;
        int get_seg_word_scan(int start, int end, uint8_t* buffer) const;
        int get_seg_word_scan_copy8(int start, int end,
                                    uint8_t* buffer) const;
        bool check_segment(int start, int end, const uint8_t* true_text,
                           uint8_t* buffer) const;
        bool check_segment_word_scan(int start, int end,
                                     const uint8_t* true_text,
                                     uint8_t* buffer) const;
        bool check_segment_word_scan_copy8(int start, int end,
                                           const uint8_t* true_text,
                                           uint8_t* buffer) const;
        uint64_t space_usage() const;
        uint64_t packed_offsets_size() const;
        uint64_t offset_low_size() const;
        uint64_t offset_presence_size() const;
        uint64_t offset_ternary_size() const;
        unsigned offset_width() const;
    protected:

    private:
        pasta::BitVector bv;
        pasta::Rank<pasta::OptimizedFor::ONE_QUERIES> *bvr;
        PackedTernaryOffsets packed_offsets;
        PackedFixedOffsets fixed_offsets;
        bool vbyte_offsets;
        const uint8_t* dictionary;
        int segments;
        int text_length;
        uint64_t common_size;
        uint64_t regular_block_part_size;
        int dictionary_blocks;
        uint64_t offset_at(int phrase) const;
        uint64_t physical_offset(uint64_t logical_offset,
                                 int phrase_begin) const;

};

#endif // RANK_HPP
