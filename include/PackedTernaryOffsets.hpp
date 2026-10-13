#pragma once

#include "pasta/bit_vector/bit_vector.hpp"
#include "pasta/bit_vector/support/rank.hpp"
#include <array>
#include <cstdint>
#include <stdexcept>
#include <vector>

// Offset representation:
//   * the low k-2 bits are stored for every offset;
//   * a bitvector marks offsets whose two high bits are nonzero;
//   * the three nonzero high values 01, 10, 11 are ternary-coded, five per byte.
class PackedTernaryOffsets {
public:
    PackedTernaryOffsets() = default;
    ~PackedTernaryOffsets() { delete high_rank_; }

    PackedTernaryOffsets(const PackedTernaryOffsets&) = delete;
    PackedTernaryOffsets& operator=(const PackedTernaryOffsets&) = delete;

    void build(const uint32_t* offsets, int count, uint64_t dictionary_size) {
        delete high_rank_;
        high_rank_ = nullptr;
        full_bits_ = 1;
        uint64_t max_offset = dictionary_size > 0 ? dictionary_size - 1 : 0;
        while ((max_offset >>= 1) != 0)
            ++full_bits_;
        low_bits_ = full_bits_ > 2 ? full_bits_ - 2 : 0;

        const uint64_t low_bit_count = static_cast<uint64_t>(count) * low_bits_;
        low_.assign(static_cast<size_t>((low_bit_count + 63) / 64), 0);
        high_present_.resize(static_cast<size_t>(count), 0);
        packed_trits_.clear();

        const uint64_t low_mask = low_bits_ == 0
            ? 0 : ((uint64_t{1} << low_bits_) - 1);
        static constexpr uint8_t powers[5] = {1, 3, 9, 27, 81};
        uint64_t nonzero_high = 0;

        for (int phrase = 0; phrase < count; ++phrase) {
            const uint64_t value = offsets[phrase];
            if (value >= dictionary_size)
                throw std::runtime_error("Phrase offset exceeds dictionary size");

            if (low_bits_ != 0) {
                const uint64_t bit_position =
                    static_cast<uint64_t>(phrase) * low_bits_;
                const size_t word = static_cast<size_t>(bit_position >> 6);
                const unsigned shift = static_cast<unsigned>(bit_position & 63);
                const uint64_t low_value = value & low_mask;
                low_[word] |= low_value << shift;
                if (shift + low_bits_ > 64)
                    low_[word + 1] |= low_value >> (64 - shift);
            }

            const uint8_t high = static_cast<uint8_t>(value >> low_bits_);
            if (high != 0) {
                high_present_[static_cast<size_t>(phrase)] = 1;
                if (nonzero_high % 5 == 0)
                    packed_trits_.push_back(0);
                packed_trits_.back() += static_cast<uint8_t>(
                    (high - 1) * powers[nonzero_high % 5]);
                ++nonzero_high;
            }
        }
        high_rank_ =
            new pasta::Rank<pasta::OptimizedFor::ONE_QUERIES>(high_present_);
    }

    uint64_t at(int phrase) const {
        uint64_t low_value = 0;
        if (low_bits_ != 0) {
            const uint64_t bit_position =
                static_cast<uint64_t>(phrase) * low_bits_;
            const size_t word = static_cast<size_t>(bit_position >> 6);
            const unsigned shift = static_cast<unsigned>(bit_position & 63);
            low_value = low_[word] >> shift;
            if (shift + low_bits_ > 64)
                low_value |= low_[word + 1] << (64 - shift);
            low_value &= (uint64_t{1} << low_bits_) - 1;
        }

        if (!high_present_[static_cast<size_t>(phrase)])
            return low_value;

        const uint64_t ordinal = high_rank_->rank1(
            static_cast<uint64_t>(phrase) + 1) - 1;
        const uint16_t decoded = trit_lookup_[packed_trits_[ordinal / 5]];
        const uint8_t high = static_cast<uint8_t>(
            ((decoded >> (2 * (ordinal % 5))) & 3U) + 1U);
        return low_value | (static_cast<uint64_t>(high) << low_bits_);
    }

    uint64_t space_usage() const {
        return static_cast<uint64_t>(low_.size()) * sizeof(uint64_t) +
               static_cast<uint64_t>(packed_trits_.size()) +
               static_cast<uint64_t>(high_present_.space_usage()) +
               (high_rank_ ? static_cast<uint64_t>(high_rank_->space_usage()) : 0);
    }

    uint64_t low_size() const {
        return static_cast<uint64_t>(low_.size()) * sizeof(uint64_t);
    }
    uint64_t presence_size() const {
        return static_cast<uint64_t>(high_present_.space_usage()) +
               (high_rank_ ? static_cast<uint64_t>(high_rank_->space_usage()) : 0);
    }
    uint64_t ternary_size() const { return packed_trits_.size(); }
    unsigned width() const { return full_bits_; }
    unsigned low_width() const { return low_bits_; }

private:
    inline static constexpr std::array<uint16_t, 243> trit_lookup_ = [] {
        std::array<uint16_t, 243> table{};
        for (unsigned value = 0; value < table.size(); ++value) {
            unsigned x = value;
            uint16_t decoded = 0;
            for (unsigned position = 0; position < 5; ++position) {
                decoded |= static_cast<uint16_t>(x % 3) << (2 * position);
                x /= 3;
            }
            table[value] = decoded;
        }
        return table;
    }();
    std::vector<uint64_t> low_;
    pasta::BitVector high_present_;
    pasta::Rank<pasta::OptimizedFor::ONE_QUERIES>* high_rank_ = nullptr;
    std::vector<uint8_t> packed_trits_;
    unsigned full_bits_ = 1;
    unsigned low_bits_ = 0;
};
