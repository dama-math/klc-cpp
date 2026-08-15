#pragma once

#include "Permutation.hpp"

/**
 * @brief Namespace for operations related to the Bruhat order on permutations.
 */
namespace Bruhat {

/**
 * @brief Reduces the permutation pair (x, w) in place.
 * 
 * This reduction applies standard transformations to simplify the pair 
 * without changing their Bruhat order relationship or associated KL polynomials.
 * Corresponds to olp_reduce_inplace in the original C algorithms.
 * 
 * @param x The first permutation (will be modified in place).
 * @param w The second permutation (will be modified in place).
 */
bool reduce(Permutation& x, Permutation& w);

/**
 * @brief Checks if the permutations (x, w) form an extremal pair.
 * 
 * @param x The first permutation.
 * @param w The second permutation.
 * @return True if (x, w) is an extremal pair, false otherwise.
 */
inline bool is_extremal_pair(const Permutation& x, const Permutation& w) {
    int n = x.size();
    if (n != w.size()) return false;

    const auto* x_ptr = x.begin();
    const auto* w_ptr = w.begin();

    // Check direct condition
    for (int i = 0; i < n - 1; ++i) {
        if (x_ptr[i] < x_ptr[i+1] && w_ptr[i] > w_ptr[i+1]) {
            return false;
        }
    }

    // Prepare arrays for the inverse check
    int8_t x_inv[16];
    int8_t w_inv[16];
    
    for (int i = 0; i < n; ++i) {
        x_inv[static_cast<int>(x_ptr[i])] = static_cast<int8_t>(i);
        w_inv[static_cast<int>(w_ptr[i])] = static_cast<int8_t>(i);
    }

    // Check inverse condition
    for (int i = 0; i < n - 1; ++i) {
        if (x_inv[i] < x_inv[i+1] && w_inv[i] > w_inv[i+1]) {
            return false;
        }
    }

    return true;
}

/**
 * @brief Checks if permutation x is less than or equal to permutation w in the Bruhat order.
 * 
 * @param x The first permutation.
 * @param w The second permutation.
 * @return True if x <= w in the Bruhat order, false otherwise.
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
    
    int8_t diff_row[16] = {0};

    const auto* x_ptr = x.begin();
    const auto* w_ptr = w.begin();

    for (int i = 0; i < n; ++i) {
        int xi = x_ptr[i];
        int wi = w_ptr[i];
        
        if (xi < wi) {
            for (int j = xi; j < wi; ++j) {
                diff_row[j]++;
            }
        } else if (xi > wi) {
            for (int j = wi; j < xi; ++j) {
                if (--diff_row[j] < 0) return false;
            }
        }
    }

    return true;
}

} // namespace Bruhat
