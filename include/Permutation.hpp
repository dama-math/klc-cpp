#pragma once

#include <string>
#include <vector>
#include <cstdint>
#include <array>
#include <stdexcept>
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
    Permutation() : n(0), length_(-1) {}

    /** @brief Construct directly from raw pointer (optimized). 
     *  Assume data checks MAX_SIZE externally or trusted source. 
     */
    Permutation(const ValueType* data_ptr, int n);

    /** @brief Construct from int pointer (requires cast/copy). */
    Permutation(const int* data_ptr, int n);

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

    /** @retutn Constant iterator to the beginning. */
    const ValueType* begin() const { return data.data(); }

    /** @retutn Constant iterator to the end. */
    const ValueType* end() const { return data.data() + n; }

    /** @return The size of the permutation. */
    int size() const { return n; }

    /** @return The Coxeter length of the permutation. */
    int length() const;

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

    /** @brief Updates the cached size based on internal n (no-op now, n is source of truth). */
    void update_size() { /* n is maintained directly */ }

    /** @brief Resets the cached length value. */
    void reset_length_cache() { length_ = -1; }
    
    /** @brief Sets size manually (use with caution when writing to data_ptr_mut directly). */
    void set_size(int new_n) { 
        if (new_n > MAX_SIZE) throw std::length_error("Permutation size exceeds MAX_SIZE");
        n = new_n; 
    }

    /** @return A packed representation of the permutation as a PackedInt. */
    PackedInt pack() const;

    /**
     * @brief Unpacks a permutation from a packed integer.
     * @param packed Packed representation.
     * @param n Expected size of the permutation.
     * @return The unpacked vector (helper for construction).
     */
    static std::vector<int> unpack_to_vector(PackedInt packed, int n);
    
    /** @brief Factory from packed. */
    static Permutation from_packed(PackedInt packed, int n);

    /** @return The inverse permutation. */
    Permutation inverse() const;

    /**
     * @brief Multiplies the permutation by a simple transposition (i, i+1).
     * @param i Index of the transposition.
     * @return The resulting permutation.
     */
    Permutation multiply_simple_transposition(int i) const;

    /** @brief Equality operator. */
    bool operator==(const Permutation& rhs) const;
    
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
