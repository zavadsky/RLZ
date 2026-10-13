#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

class PackedFixedOffsets {
public:
    void build(const uint32_t* offsets, int count, uint64_t dictionary_size) {
        width_ = 1;
        uint64_t maximum = dictionary_size > 0 ? dictionary_size - 1 : 0;
        while ((maximum >>= 1) != 0)
            ++width_;
        const uint64_t bits = static_cast<uint64_t>(count) * width_;
        words_.assign(static_cast<size_t>((bits + 63) / 64), 0);
        for (int phrase = 0; phrase < count; ++phrase) {
            const uint64_t value = offsets[phrase];
            if (value >= dictionary_size)
                throw std::runtime_error("Phrase offset exceeds dictionary size");
            const uint64_t bit = static_cast<uint64_t>(phrase) * width_;
            const size_t word = static_cast<size_t>(bit >> 6);
            const unsigned shift = static_cast<unsigned>(bit & 63);
            words_[word] |= value << shift;
            if (shift + width_ > 64)
                words_[word + 1] |= value >> (64 - shift);
        }
    }

    uint64_t at(int phrase) const {
        const uint64_t bit = static_cast<uint64_t>(phrase) * width_;
        const size_t word = static_cast<size_t>(bit >> 6);
        const unsigned shift = static_cast<unsigned>(bit & 63);
        uint64_t value = words_[word] >> shift;
        if (shift + width_ > 64)
            value |= words_[word + 1] << (64 - shift);
        const uint64_t mask = width_ == 64 ? ~uint64_t{0}
                                           : ((uint64_t{1} << width_) - 1);
        return value & mask;
    }

    uint64_t space_usage() const { return words_.size() * sizeof(uint64_t); }
    unsigned width() const { return width_; }

private:
    std::vector<uint64_t> words_;
    unsigned width_ = 1;
};
