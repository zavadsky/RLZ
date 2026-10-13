#ifndef RLZ_RANK_SDSL_HPP
#define RLZ_RANK_SDSL_HPP
#include <cstdint>
#include <vector>
#include "PackedTernaryOffsets.hpp"
#include "PackedFixedOffsets.hpp"


class RLZ_Rank_ef
{
    public:
        RLZ_Rank_ef(const uint32_t* lengths, const uint32_t* offsets,
                    const uint8_t* dictionary, int segments,
                    int text_length, uint64_t dictionary_size,
                    uint64_t common_size = 0,
                    uint64_t regular_block_part_size = 0,
                    int dictionary_blocks = 1,
                    bool vbyte_offsets = true);
        virtual ~RLZ_Rank_ef();
        int get_seg(int start, int end, uint8_t* buffer) const;
        bool check_segment(int start, int end, const uint8_t* true_text,
                           uint8_t* buffer) const;
        uint64_t space_usage() const;
        uint64_t packed_offsets_size() const;
        unsigned offset_width() const;

    private:
        struct EFData;
        EFData* ef_data;
        PackedTernaryOffsets packed_offsets;
        PackedFixedOffsets fixed_offsets;
        bool vbyte_offsets;
        const uint8_t* dictionary;
        int segments;
        int text_length;
        uint64_t common_size, regular_block_part_size;
        int dictionary_blocks;
        uint64_t offset_at(int phrase) const;
        uint64_t physical_offset(uint64_t logical_offset,
                                 int phrase_begin) const;
        uint64_t low_at(int phrase) const;
        uint64_t value_at(int phrase, uint64_t high_position) const;
        uint64_t next_high_one(uint64_t position) const;
        int predecessor(int position, uint64_t& high_position) const;
};



#endif // RLZ_RANK_SDSL_HPP
