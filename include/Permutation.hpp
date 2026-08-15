#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <array>
#include <stdexcept>
#include <bit>
#include <algorithm>
#include "Config.hpp"

/**
 * @brief Represents a permutation.
 * 
 * Optimized to use a fixed-size array for storage (up to size 32).
 * This makes the class a trivial type (POD-like) for efficient copy and passing.
 */
class Permutation {
public:
    static constexpr int MAX_SIZE = MAX_PACKED_SIZE;
    using ValueType = int8_t; // 8-bit integer is enough for n <= 24

    /** @brief Default constructor (empty permutation). */
    Permutation() : n(0) {}

    /** @brief Construct directly from raw pointer (optimized). 
     *  Assume data checks MAX_SIZE externally or trusted source. 
     */
    Permutation(const ValueType* data_ptr, int n_val) {
        n = n_val;
        std::copy(data_ptr, data_ptr + n, data.begin());
    }

    /** @brief Construct from int pointer (requires cast/copy). */
    Permutation(const int* data_ptr, int n_val) {
        n = n_val;
        for (int i = 0; i < n; ++i) {
            data[i] = static_cast<ValueType>(data_ptr[i]);
        }
    }

    /**
     * @brief Constructs a permutation from a string (one-line notation).
     * @param one_line_notation String representing the permutation (e.g., "123").
     */
    Permutation(const std::string& one_line_notation);

    /**
     * @brief Constructs a permutation from a vector of integers.
     * @param data_vec Vector containing the permutation elements.
     */
    Permutation(const std::vector<int>& data_vec);

    /** @return Constant iterator to the beginning. */
    const ValueType* begin() const { return data.data(); }

    /** @return Constant iterator to the end. */
    const ValueType* end() const { return data.data() + n; }

    /** @return The size of the permutation. */
    int size() const { return n; }

    /** @return The Coxeter length of the permutation. */
    int length() const {
        if (length_ != -1) return length_;
        int inversions = 0;
        uint32_t seen = 0;
        for (int i = 0; i < n; ++i) {
            int v = data[i];
            if (v < 31) {
                inversions += std::popcount(seen >> (v + 1));
            }
            seen |= (1u << v);
        }
        length_ = inversions;
        return length_;
    }

    /** 
     * @brief Direct access to elements.
     * @param i Index.
     */
    int operator[](int i) const { return data[i]; }

    /**
     * @brief Returns a copy of data as std::vector<int>.
     * Primarily for legacy compatibility and printing.
     */
    std::vector<int> to_vector() const;
    
    // Low-level access for optimization (e.g. reductions)
    ValueType* data_ptr_mut() { return data.data(); }
    const ValueType* data_ptr() const { return data.data(); }

    /** @brief Resets the cached length value. */
    void reset_length_cache() { length_ = -1; }
    
    bool has_cached_length() const { return length_ != -1; }
    
    void set_cached_length(int len) { length_ = len; }

    /** @brief Sets size manually (use with caution when writing to data_ptr_mut directly). */
    void set_size(int new_n) { 
        if (new_n > MAX_SIZE) throw std::length_error("Permutation size exceeds MAX_SIZE");
        n = new_n; 
    }

    /** @return A packed representation of the permutation as a PackedInt. */
    PackedInt pack() const {
        PackedInt packed = 0;
        for (int i = 0; i < n; ++i) {
            packed |= ((PackedInt)(data[i] & PACKED_MASK) << (BITS_PER_ELEMENT * i));
        }
        return packed;
    }

    /** @return The inverse permutation. */
    Permutation inverse() const {
        Permutation new_perm;
        new_perm.n = n;
        for (int i = 0; i < n; ++i) {
            new_perm.data[static_cast<int>(data[i])] = static_cast<ValueType>(i);
        }
        if (this->length_ != -1) {
            new_perm.length_ = this->length_;
        }
        return new_perm;
    }

    /** @return The Dynkin diagram automorphism permutation theta(w) = w_0 w w_0. */
    Permutation theta() const {
        Permutation new_perm;
        new_perm.n = n;
        for (int i = 0; i < n; ++i) {
            new_perm.data[i] = static_cast<ValueType>((n - 1) - data[n - 1 - i]);
        }
        if (this->length_ != -1) {
            new_perm.length_ = this->length_;
        }
        return new_perm;
    }

    /**
     * @brief Multiplies the permutation by a simple transposition (i, i+1).
     * @param i Index of the transposition.
     * @return The resulting permutation.
     */
    Permutation multiply_simple_transposition(int i) const {
        Permutation new_perm = *this;
        std::swap(new_perm.data[i], new_perm.data[i + 1]);
        if (this->length_ != -1) {
            new_perm.length_ = (data[i] < data[i + 1]) ? (this->length_ + 1) : (this->length_ - 1);
        }
        return new_perm;
    }

    /** @brief Equality operator. */
    bool operator==(const Permutation& rhs) const {
        if (n != rhs.n) return false;
        for (int i = 0; i < n; ++i) {
            if (data[i] != rhs.data[i]) return false;
        }
        return true;
    }
    
    /** @return String representation in one-line notation. */
    std::string to_string() const;

    /**
     * @brief Converts a data vector to a one-line notation string.
     * @param data Data vector.
     * @return One-line notation string.
     */
    static std::string to_one_line_notation(const std::vector<int>& data);

private:
    std::array<ValueType, MAX_SIZE> data{}; ///< 0-indexed data.
    int n;                 ///< Size of the permutation.
    mutable int length_ = -1; ///< Cached length.

    /**
     * @brief Parses a one-line notation string.
     * @param notation One-line notation string.
     * @return Data vector.
     */
    static std::vector<int> from_one_line_notation(const std::string& notation);
};

