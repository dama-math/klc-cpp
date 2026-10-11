#pragma once

#include "Permutation.hpp"

/**
 * @brief Algorithms for the strong Bruhat order and Warrington's pair reductions.
 */
namespace Bruhat {

/**
 * @brief Reduces the permutation pair (x, w) in place.
 * 
 * Applies Warrington's extremal-pair swaps and tilde-reduction (eliminating
 * fixed values with zero interval overlap). Preserves the Kazhdan-Lusztig polynomial
 * P_{x, w}(q) while simplifying the permutations.
 * 
 * @return True if (x, w) was originally an extremal pair before reduction.
 */
bool reduce(Permutation& x, Permutation& w);

/**
 * @brief Checks whether (x, w) forms an extremal pair.
 * 
 * An extremal pair satisfies:
 * 1. For all i: if w has a right descent at i (w[i] > w[i+1]), x also has a descent (x[i] > x[i+1]).
 * 2. For all i: if w has a left descent at i (w^{-1}[i] > w^{-1}[i+1]), x also has a left descent.
 */
inline bool is_extremal_pair(const Permutation& x, const Permutation& w) {
    int n = x.size();
    if (n != w.size()) return false;

    const auto* x_ptr = x.begin();
    const auto* w_ptr = w.begin();

    // Check right descent condition: x must have a descent whenever w has a descent
    for (int i = 0; i < n - 1; ++i) {
        if (x_ptr[i] < x_ptr[i + 1] && w_ptr[i] > w_ptr[i + 1]) {
            return false;
        }
    }

    // Prepare inverses for left descent condition check
    int8_t x_inv[Permutation::MAX_SIZE];
    int8_t w_inv[Permutation::MAX_SIZE];
    for (int i = 0; i < n; ++i) {
        x_inv[static_cast<int>(x_ptr[i])] = static_cast<int8_t>(i);
        w_inv[static_cast<int>(w_ptr[i])] = static_cast<int8_t>(i);
    }

    for (int i = 0; i < n - 1; ++i) {
        if (x_inv[i] < x_inv[i + 1] && w_inv[i] > w_inv[i + 1]) {
            return false;
        }
    }

    return true;
}

struct DiffIntervalTable {
    alignas(16) int8_t data[16][16][16];
    constexpr DiffIntervalTable() : data{} {
        for (int a = 0; a < 16; ++a) {
            for (int b = 0; b < 16; ++b) {
                for (int k = 0; k < 16; ++k) {
                    data[a][b][k] = (k >= a && k < b) ? 1 : 0;
                }
            }
        }
    }
};

inline constexpr DiffIntervalTable DIFF_INTERVAL_TABLE_HOLDER;
#define DIFF_INTERVAL_TABLE DIFF_INTERVAL_TABLE_HOLDER.data

/**
 * @brief Checks if permutation x is less than or equal to w in the strong Bruhat order (x <= w).
 * 
 * Uses SIMD parallel interval additions/subtractions on 16-byte diff vector.
 */
inline bool is_bruhat_less_than(const Permutation& x, const Permutation& w) {
    if (x.size() != w.size()) return false;
    int n = x.size();
    if (n == 0) return true;
    
    int x_len = x.length();
    if (x_len == 0) return true;
    int w_len = w.length();
    if (x_len > w_len) return false;
    if (x_len == w_len) return x == w;
    
    const auto* x_ptr = x.begin();
    const auto* w_ptr = w.begin();

    __m128i diff = _mm_setzero_si128();

    for (int i = 0; i < n; ++i) {
        int xi = x_ptr[i];
        int wi = w_ptr[i];
        
        if (xi < wi) {
            __m128i mask = _mm_loadu_si128(reinterpret_cast<const __m128i*>(DIFF_INTERVAL_TABLE[xi][wi]));
            diff = _mm_add_epi8(diff, mask);
        } else if (xi > wi) {
            __m128i mask = _mm_loadu_si128(reinterpret_cast<const __m128i*>(DIFF_INTERVAL_TABLE[wi][xi]));
            diff = _mm_sub_epi8(diff, mask);
            if (_mm_movemask_epi8(diff) != 0) return false;
        }
    }

    return true;
}

} // namespace Bruhat
