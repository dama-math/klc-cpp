#include "Bruhat.hpp"
#include <algorithm>
#include <bit>
#include <cstdint>

namespace Bruhat {

namespace {

using ValueType = Permutation::ValueType;
constexpr ValueType REMOVED_MARKER = -1;

/**
 * @brief Enforces extremal-pair conditions on (x, w) by performing descent-preserving swaps on x.
 * 
 * If w has a descent at index i (or w^{-1} at index i) where x has an ascent,
 * swapping the corresponding elements in x preserves the KL polynomial while increasing l(x).
 * 
 * @return Total number of swaps performed on x.
 */
int make_extremal_pair_inplace(ValueType* x_data, const ValueType* w_data, int n) {
    if (n <= 1) return 0;
    
    ValueType x_inv[Permutation::MAX_SIZE];
    ValueType w_inv[Permutation::MAX_SIZE];
    
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
        
        // Check right descents
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
        
        // Check left descents (descents in inverse permutation)
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
 * @brief Identifies fixed points x[i] == w[i] that do not cross any active intervals (tilde-reduction).
 * 
 * If x[i] == w[i] and the net interval crossing count dxw is 0, element i does not affect
 * the KL polynomial and can be removed from both permutations.
 * 
 * @return Number of elements marked as REMOVED_MARKER.
 */
int tilde_reduction_inplace(ValueType* x_data, ValueType* w_data, int n) {
    if (n <= 1) return 0;
    int removed_count = 0;

    for (int i = 0; i < n; ++i) {
        if (x_data[i] == w_data[i]) {
            int x_val = x_data[i];
            int dxw = 0;
            for (int j = 0; j < i; ++j) {
                int xj = x_data[j];
                int wj = w_data[j];
                if (xj <= x_val && x_val < wj) {
                    dxw++;
                } else if (wj <= x_val && x_val < xj) {
                    dxw--;
                }
            }
            if (dxw == 0) {
                x_data[i] = REMOVED_MARKER;
                w_data[i] = REMOVED_MARKER;
                removed_count++;
            }
        }
    }
    return removed_count;
}

/**
 * @brief Compacts the permutation arrays after removing elements, re-indexing values to 0..n'-1.
 */
void flatten_inplace(ValueType* x_data, ValueType* w_data, int& n) {
    if (n == 0) return;

    uint32_t active_vals = 0;
    for (int i = 0; i < n; ++i) {
        if (x_data[i] != REMOVED_MARKER) {
            active_vals |= (1u << x_data[i]);
            active_vals |= (1u << w_data[i]);
        }
    }
    
    ValueType rank_table[Permutation::MAX_SIZE];
    int r = 0;
    uint32_t m = active_vals;
    while (m) {
        int val = std::countr_zero(m);
        rank_table[val] = static_cast<ValueType>(r++);
        m &= m - 1;
    }

    int write_idx = 0;
    for (int i = 0; i < n; ++i) {
        if (x_data[i] != REMOVED_MARKER) {
            x_data[write_idx] = rank_table[static_cast<int>(x_data[i])];
            w_data[write_idx] = rank_table[static_cast<int>(w_data[i])];
            ++write_idx;
        }
    }
    n = write_idx;
}

} // anonymous namespace

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
        int swaps = make_extremal_pair_inplace(x_data, w_data, n);
        if (first_iteration && swaps > 0) {
            originally_ep = false;
        }
        first_iteration = false;
        total_swaps += swaps;

        int removed = tilde_reduction_inplace(x_data, w_data, n);
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