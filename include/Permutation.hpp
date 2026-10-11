#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <array>
#include <stdexcept>
#include <bit>
#include <cstring>
#include <immintrin.h>
#include "Config.hpp"

/**
 * @brief Represents a permutation in the symmetric group S_n.
 * 
 * Stored as a 0-indexed flat array of up to MAX_SIZE (16 in 64-bit mode, 24 in 128-bit mode)
 * elements. Designed as a lightweight, trivially copyable value type.
 */
class Permutation {
public:
    static constexpr int MAX_SIZE = MAX_PACKED_SIZE;
    using ValueType = int8_t;

    /** @brief Constructs an empty permutation of size 0. */
    Permutation() : n(0) {}

    /** @brief Constructs a permutation directly from raw array buffer. */
    Permutation(const ValueType* data_ptr, int n_val) {
        if (n_val > MAX_SIZE) {
            throw std::length_error("Permutation size exceeds MAX_SIZE");
        }
        n = n_val;
        std::copy(data_ptr, data_ptr + n, data.begin());
    }

    /** @brief Constructs a permutation from a vector of integers. */
    Permutation(const std::vector<int>& data_vec);

    /** @brief Constructs a permutation from a 0-indexed one-line notation string (e.g. "0213"). */
    Permutation(const std::string& one_line_notation);

    const ValueType* begin() const { return data.data(); }
    const ValueType* end() const { return data.data() + n; }

    int size() const { return n; }

    /**
     * @brief Returns the Coxeter length (number of inversions) of the permutation.
     * Computes inversions using hardware bit operations if not cached.
     */
    int length() const {
        if (length_ != -1) return length_;
        int inversions = 0;
        uint32_t seen = 0;
        for (int i = 0; i < n; ++i) {
            int v = data[i];
            inversions += std::popcount(seen >> (v + 1));
            seen |= (1u << v);
        }
        length_ = inversions;
        return length_;
    }

    int operator[](int i) const { return data[i]; }

    /** @brief Direct buffer access for in-place algorithms (e.g., Bruhat reduction). */
    ValueType* data_ptr_mut() { return data.data(); }
    const ValueType* data_ptr() const { return data.data(); }

    void reset_length_cache() { length_ = -1; }
    bool has_cached_length() const { return length_ != -1; }
    void set_cached_length(int len) { length_ = len; }

    void set_size(int new_n) {
        if (new_n > MAX_SIZE) {
            throw std::length_error("Permutation size exceeds MAX_SIZE");
        }
        n = new_n;
    }

    /**
     * @brief Packs permutation elements into an integer for fast hash table lookups.
     */
    PackedInt pack() const {
        PackedInt packed = 0;
        for (int i = 0; i < n; ++i) {
            packed |= ((PackedInt)(data[i] & PACKED_MASK) << (BITS_PER_ELEMENT * i));
        }
        return packed;
    }

    /** @brief Returns the group-theoretic inverse permutation w^{-1}. */
    Permutation inverse() const {
        Permutation inv;
        inv.n = n;
        for (int i = 0; i < n; ++i) {
            inv.data[static_cast<int>(data[i])] = static_cast<ValueType>(i);
        }
        inv.length_ = this->length_;
        return inv;
    }

    /**
     * @brief Returns the permutation multiplied on the right by simple transposition s_i = (i, i+1).
     * Efficiently updates the Coxeter length: l(w * s_i) = l(w) +- 1.
     */
    Permutation multiply_simple_transposition(int i) const {
        Permutation result = *this;
        std::swap(result.data[i], result.data[i + 1]);
        if (this->length_ != -1) {
            result.length_ = (data[i] < data[i + 1]) ? (this->length_ + 1) : (this->length_ - 1);
        }
        return result;
    }

    bool operator==(const Permutation& rhs) const {
        if (n != rhs.n) return false;
        for (int i = 0; i < n; ++i) {
            if (data[i] != rhs.data[i]) return false;
        }
        return true;
    }

    /** @brief Converts the permutation to one-line notation string. */
    std::string to_string() const;

    /** @brief Converts elements to std::vector<int> for convenience. */
    std::vector<int> to_vector() const;

private:
    std::array<ValueType, MAX_SIZE> data{};
    int n = 0;
    mutable int length_ = -1;

    void parse_one_line_notation(const std::string& notation);
};
