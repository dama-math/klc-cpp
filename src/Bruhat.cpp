#include "Bruhat.hpp"
#include <vector>
#include <numeric>
#include <algorithm>
#include <iostream>

namespace Bruhat {

// Anonymous namespace for internal helper functions and logic
namespace {

using ValueType = int8_t; // Must match Permutation::ValueType
constexpr ValueType REMOVED = -1;

/**
 * @brief Attempts to apply an extremal pair reduction in place.
 * 
 * Modifies the data arrays directly to transform the pair into an extremal pair
 * if possible, using inversion swaps.
 * 
 * @param x_data Array pointer to the inner permutation data.
 * @param w_data Array pointer to the outer permutation data.
 * @param n Size of the permutations.
 * @return True if the inner permutation 'x' was modified, false otherwise.
 */
int make_ep_inplace(ValueType* x_data, ValueType* w_data, int n) {
    if (n <= 1) return 0;
    
    ValueType x_inv[32];
    ValueType w_inv[32];
    
    for (int i = 0; i < n; ++i) {
        x_inv[static_cast<int>(x_data[i])] = static_cast<ValueType>(i);
        w_inv[static_cast<int>(w_data[i])] = static_cast<ValueType>(i);
    }

    uint32_t w_desc = 0;
    uint32_t w_inv_desc = 0;
    for (int i = 0; i < n - 1; ++i) {
        if (w_data[i] > w_data[i + 1]) w_desc |= (1u << i);
        if (w_inv[i] > w_inv[i + 1]) w_inv_desc |= (1u << i);
    }

    if (w_desc == 0 && w_inv_desc == 0) return 0;
    
    bool changed = true;
    int swaps_made = 0;
    
    while (changed) {
        changed = false;
        
        uint32_t mask = w_desc;
        while (mask) {
            int i = std::countr_zero(mask);
            mask &= mask - 1;
            if (x_data[i] < x_data[i + 1]) {
                int val1 = x_data[i];
                int val2 = x_data[i + 1];
                std::swap(x_data[i], x_data[i + 1]);
                std::swap(x_inv[val1], x_inv[val2]);
                changed = true;
                swaps_made++;
            }
        }
        
        uint32_t inv_mask = w_inv_desc;
        while (inv_mask) {
            int i = std::countr_zero(inv_mask);
            inv_mask &= inv_mask - 1;
            if (x_inv[i] < x_inv[i + 1]) {
                std::swap(x_data[x_inv[i]], x_data[x_inv[i + 1]]);
                std::swap(x_inv[i], x_inv[i + 1]);
                changed = true;
                swaps_made++;
            }
        }
    }
    return swaps_made;
}

/**
 * @brief Analyzes and marks elements that can be eliminated (tilde reduction).
 * 
 * Elements marked for elimination are replaced with -1.
 * 
 * @param x_data Array pointer to the inner permutation data.
 * @param w_data Array pointer to the outer permutation data.
 * @param n Size of the permutations.
 * @return The number of elements successfully marked for removal.
 */
int tilde_inplace(ValueType* x_data, ValueType* w_data, int n) {
    if (n <= 1) return 0;
    int removed_count = 0;

    for (int i = 0; i < n; ++i) {
        if (x_data[i] == w_data[i]) {
            int x = x_data[i];
            int dxw = 0;
            // Compute dxw_row[x] by looking at all previous elements
            for (int j = 0; j < i; ++j) {
                int xj = x_data[j];
                int wj = w_data[j];
                if (xj <= x && x < wj) {
                    dxw++;
                } else if (wj <= x && x < xj) {
                    dxw--;
                }
            }
            // Elements that are identical and have no overlapping intervals can be removed
            if (dxw == 0) {
                x_data[i] = REMOVED;
                w_data[i] = REMOVED;
                removed_count++;
            }
        }
    }
    return removed_count;
}

/**
 * @brief Flattens the permutations after tilde reduction.
 * 
 * Compacts the data arrays by removing elements marked as -1 and re-indexing
 * the remaining elements sequentially.
 * 
 * @param x_data Array pointer to the inner permutation data.
 * @param w_data Array pointer to the outer permutation data.
 * @param n Reference to the permutation size (updated upon completion).
 */
void flatten_inplace(ValueType* x_data, ValueType* w_data, int& n) {
    if (n == 0) return;

    uint32_t active_vals = 0;
    for (int i = 0; i < n; ++i) {
        if (x_data[i] != REMOVED) {
            active_vals |= (1u << x_data[i]);
            active_vals |= (1u << w_data[i]);
        }
    }
    
    ValueType rank_table[32];
    int r = 0;
    uint32_t m = active_vals;
    while (m) {
        int val = std::countr_zero(m);
        rank_table[val] = static_cast<ValueType>(r++);
        m &= m - 1;
    }

    int write_idx = 0;
    for (int i = 0; i < n; ++i) {
        if (x_data[i] != REMOVED) {
            x_data[write_idx] = rank_table[static_cast<int>(x_data[i])];
            w_data[write_idx] = rank_table[static_cast<int>(w_data[i])];
            ++write_idx;
        }
    }
    n = write_idx; // Update physical size
}

} // namespace

// --- Public API Functions ---

bool reduce(Permutation& x, Permutation& w) {
    auto* x_data = x.data_ptr_mut();
    auto* w_data = w.data_ptr_mut();
    int n = x.size();
    if (n <= 1) return true;

    bool originally_ep = true;
    bool first_iteration = true;
    int total_swaps = 0;
    bool any_removals = false;
    int orig_x_len = x.has_cached_length() ? x.length() : -1;

    while (true) {
        int swaps = make_ep_inplace(x_data, w_data, n);
        if (first_iteration && swaps > 0) {
            originally_ep = false;
        }
        first_iteration = false;
        total_swaps += swaps;

        int removed = tilde_inplace(x_data, w_data, n);
        if (removed > 0) {
            any_removals = true;
            flatten_inplace(x_data, w_data, n);
            x.set_size(n);
            w.set_size(n);
        }

        if (swaps == 0 && removed == 0) {
            break;
        }
    }

    if (any_removals) {
        x.reset_length_cache();
        w.reset_length_cache();
    } else if (total_swaps > 0) {
        if (orig_x_len != -1) {
            x.set_cached_length(orig_x_len + total_swaps);
        } else {
            x.reset_length_cache();
        }
    }

    return originally_ep;
}

} // namespace Bruhat