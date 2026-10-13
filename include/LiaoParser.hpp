#pragma once

#include <algorithm>
#include <array>
#include <bit>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <vector>

#include "libsais64.h"
#include <sdsl/rmq_support.hpp>

// Greedy RLZ parsing kernel adapted from mpetri/rlz-store's
// dict_index_csa.hpp and factor_select_first. As in its fast flat FM-index
// design, each dictionary symbol has its own BWT bitvector and rank directory.
// Index construction, allocation, and verification are separate from parse().
class LiaoParser {
public:
    explicit LiaoParser(std::size_t text_size)
        : offsets_(text_size), lengths_(text_size) {}

    void build_index(const uint8_t* dictionary,
                     std::size_t dictionary_size,
                     bool build_leftmost_rmq = true) {
        dictionary_size_ = dictionary_size;
        if (!dictionary || dictionary_size == 0)
            throw std::invalid_argument("Liao parser requires a dictionary");
        if (std::find(dictionary, dictionary + dictionary_size, uint8_t{0}) !=
            dictionary + dictionary_size)
            throw std::invalid_argument(
                "The rlz-store FM index reserves byte 0 as its sentinel");

        std::vector<uint8_t> reversed(dictionary_size + 1, 0);
        std::reverse_copy(dictionary, dictionary + dictionary_size,
                          reversed.begin());
        suffix_array_.resize(reversed.size());
        if (libsais64(reversed.data(), suffix_array_.data(),
                      static_cast<int64_t>(reversed.size()), 0, nullptr) < 0)
            throw std::runtime_error("Cannot construct Liao suffix array");

        // A phrase occurrence represented by suffix-array position i starts
        // at dictionary_size_ - (suffix_array_[i] + phrase_length) in the
        // original (non-reversed) dictionary.  Consequently, the leftmost
        // occurrence is obtained from the maximum suffix-array value in the
        // final FM interval.  The succinct RMQ occupies 2n+o(n) bits and
        // returns that maximum's position in O(1) time.
        has_leftmost_rmq_ = build_leftmost_rmq;
        if (has_leftmost_rmq_)
            maximum_suffix_rmq_ =
                sdsl::rmq_succinct_sct<false>(&suffix_array_);

        std::array<uint64_t, 256> counts{};
        for (uint8_t symbol : reversed)
            ++counts[symbol];
        uint64_t cumulative = 0;
        for (std::size_t symbol = 0; symbol < 256; ++symbol) {
            C_[symbol] = cumulative;
            cumulative += counts[symbol];
        }
        C_[256] = cumulative;

        alphabet_size_ = 0;
        symbol_row_.fill(-1);
        for (std::size_t symbol = 1; symbol < 256; ++symbol)
            if (counts[symbol] != 0)
                symbol_row_[symbol] = static_cast<int>(alphabet_size_++);

        words_ = (reversed.size() + 63) / 64;
        bits_.assign(alphabet_size_ * words_, 0);
        for (std::size_t bwt_position = 0;
             bwt_position < suffix_array_.size(); ++bwt_position) {
            const int64_t suffix = suffix_array_[bwt_position];
            const uint8_t symbol = suffix == 0
                ? uint8_t{0}
                : reversed[static_cast<std::size_t>(suffix - 1)];
            const int row = symbol_row_[symbol];
            if (row >= 0)
                bits_[static_cast<std::size_t>(row) * words_ +
                      bwt_position / 64] |= uint64_t{1} << (bwt_position % 64);
        }

        rank_prefix_.assign(alphabet_size_ * (words_ + 1), 0);
        for (std::size_t row = 0; row < alphabet_size_; ++row) {
            uint32_t sum = 0;
            for (std::size_t word = 0; word < words_; ++word) {
                rank_prefix_[row * (words_ + 1) + word] = sum;
                sum += static_cast<uint32_t>(std::popcount(
                    bits_[row * words_ + word]));
            }
            rank_prefix_[row * (words_ + 1) + words_] = sum;
        }
    }

    std::size_t parse(const uint8_t* text, std::size_t text_size,
                      bool leftmost_offsets = true) {
        if (!text && text_size != 0)
            throw std::invalid_argument("Liao parser received a null text");
        if (leftmost_offsets && !has_leftmost_rmq_)
            throw std::logic_error(
                "Leftmost Liao parsing requires an RMQ-enabled index");
        std::size_t position = 0;
        std::size_t factors = 0;
        while (position < text_size) {
            uint64_t sp = 0;
            uint64_t ep = suffix_array_.size() - 1;
            std::size_t length = 0;
            while (position + length < text_size) {
                const uint8_t symbol = text[position + length];
                if (symbol_row_[symbol] < 0)
                    break;
                const uint64_t next_sp = C_[symbol] + rank(symbol, sp);
                const uint64_t next_ep = C_[symbol] + rank(symbol, ep + 1) - 1;
                if (next_ep < next_sp)
                    break;
                sp = next_sp;
                ep = next_ep;
                ++length;
            }
            if (length == 0)
                throw std::runtime_error(
                    "The dictionary does not contain a text byte");

            const std::size_t selected_position = leftmost_offsets
                ? static_cast<std::size_t>(maximum_suffix_rmq_(sp, ep))
                : static_cast<std::size_t>(sp);
            const uint64_t offset = dictionary_size_ -
                (static_cast<uint64_t>(suffix_array_[selected_position]) +
                 length);
            offsets_[factors] = static_cast<uint32_t>(offset);
            lengths_[factors] = static_cast<uint32_t>(length);
            ++factors;
            position += length;
        }
        factor_count_ = factors;
        return factors;
    }

    bool verify(const uint8_t* text, std::size_t text_size,
                const uint8_t* dictionary) const {
        std::size_t text_position = 0;
        for (std::size_t i = 0; i < factor_count_; ++i) {
            const std::size_t offset = offsets_[i];
            const std::size_t length = lengths_[i];
            if (offset + length > dictionary_size_ ||
                text_position + length > text_size ||
                !std::equal(dictionary + offset, dictionary + offset + length,
                            text + text_position))
                return false;
            text_position += length;
        }
        return text_position == text_size;
    }

    const uint32_t* offsets() const { return offsets_.data(); }
    const uint32_t* lengths() const { return lengths_.data(); }
    std::size_t factor_count() const { return factor_count_; }

private:
    uint64_t rank(uint8_t symbol, uint64_t end) const {
        const std::size_t row = static_cast<std::size_t>(symbol_row_[symbol]);
        const std::size_t word = static_cast<std::size_t>(end / 64);
        const unsigned remainder = static_cast<unsigned>(end % 64);
        uint64_t result = rank_prefix_[row * (words_ + 1) + word];
        if (remainder != 0 && word < words_) {
            const uint64_t mask = (uint64_t{1} << remainder) - 1;
            result += std::popcount(bits_[row * words_ + word] & mask);
        }
        return result;
    }

    std::size_t dictionary_size_ = 0;
    std::vector<int64_t> suffix_array_;
    sdsl::rmq_succinct_sct<false> maximum_suffix_rmq_;
    bool has_leftmost_rmq_ = false;
    std::array<uint64_t, 257> C_{};
    std::array<int, 256> symbol_row_{};
    std::size_t alphabet_size_ = 0;
    std::size_t words_ = 0;
    std::vector<uint64_t> bits_;
    std::vector<uint32_t> rank_prefix_;
    std::vector<uint32_t> offsets_;
    std::vector<uint32_t> lengths_;
    std::size_t factor_count_ = 0;
};
