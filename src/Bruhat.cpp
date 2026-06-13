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
    bool changed = true;
    int swaps_made = 0;
    
    ValueType x_inv[32];
    ValueType w_inv[32];
    
    for (int i = 0; i < n; ++i) {
        x_inv[x_data[i]] = i;
        w_inv[w_data[i]] = i;
    }
    
    while (changed) {
        changed = false;
        
        for (int i = 0; i < n - 1; ++i) {
            if (w_data[i] > w_data[i + 1] && x_data[i] < x_data[i + 1]) {
                int val1 = x_data[i];
                int val2 = x_data[i + 1];
                std::swap(x_data[i], x_data[i + 1]);
                std::swap(x_inv[val1], x_inv[val2]);
                changed = true;
                swaps_made++;
            }
        }
        
        for (int i = 0; i < n - 1; ++i) {
            if (w_inv[i] > w_inv[i + 1] && x_inv[i] < x_inv[i + 1]) {
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
    
    int write_idx = 0;
    for (int i = 0; i < n; ++i) {
        if (x_data[i] != REMOVED) {
            x_data[write_idx] = static_cast<ValueType>(__builtin_popcount(active_vals & ((1u << x_data[i]) - 1)));
            w_data[write_idx] = static_cast<ValueType>(__builtin_popcount(active_vals & ((1u << w_data[i]) - 1)));
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

    bool changed = true;
    bool modified = false;
    bool originally_ep = true;
    bool first_iteration = true;
    
    int total_x_swaps = 0;
    bool tilde_applied = false;
    
    int initial_x_len = x.has_cached_length() ? x.length() : -1;
    
    while (changed) {
        int swaps = make_ep_inplace(x_data, w_data, n);
        if (first_iteration && swaps > 0) {
            originally_ep = false;
        }
        first_iteration = false;
        
        if (swaps > 0) {
            modified = true;
            total_x_swaps += swaps;
        }
        
        int removed_count = tilde_inplace(x_data, w_data, n);
        if (removed_count > 0) {
            flatten_inplace(x_data, w_data, n);
            // Updating size in wrapper Permutation objects
            x.set_size(n);
            w.set_size(n);
            changed = true;
            modified = true;
            tilde_applied = true;
        } else {
            changed = false;
        }
    }
    
    // Conservatively invalidate length caches for both permutations
    // since make_ep_inplace can modify x data and tilde/flatten can modify both.
    if (modified) {
        if (tilde_applied) {
            x.reset_length_cache();
            w.reset_length_cache();
        } else {
            if (initial_x_len != -1) {
                x.set_cached_length(initial_x_len + total_x_swaps);
            } else {
                x.reset_length_cache();
            }
        }
    }
    return originally_ep;
}

bool is_extremal_pair(const Permutation& x, const Permutation& w) {
    int n = x.size();
    if (n != w.size()) return false;

    // Check direct condition using operator[]
    for (int i = 0; i < n - 1; ++i) {
        if (x[i] < x[i+1] && w[i] > w[i+1]) {
            return false;
        }
    }

    // Prepare arrays for the inverse check
    int x_inv[32];
    int w_inv[32];
    
    for (int i = 0; i < n; ++i) {
        x_inv[x[i]] = i;
        w_inv[w[i]] = i;
    }

    // Check inverse condition
    for (int i = 0; i < n - 1; ++i) {
        if (x_inv[i] < x_inv[i+1] && w_inv[i] > w_inv[i+1]) {
            return false;
        }
    }

    return true;
}

bool is_bruhat_less_than(const Permutation& x, const Permutation& w) {
    if (x.size() != w.size()) return false;
    int n = x.size();
    if (n == 0) return true;
    
    int8_t diff_row[32] = {0};

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

bool is_extremal_pair(const std::vector<int>& x_data, const std::vector<int>& w_data) {
    int n = x_data.size();
    if (n != (int)w_data.size()) return false;

    for (int i = 0; i < n - 1; ++i) {
        if (x_data[i] < x_data[i+1] && w_data[i] > w_data[i+1]) return false;
    }

    std::vector<int> x_inv(n), w_inv(n);
    for (int i = 0; i < n; ++i) {
        x_inv[x_data[i]] = i;
        w_inv[w_data[i]] = i;
    }
    
    for (int i = 0; i < n - 1; ++i) {
        if (x_inv[i] < x_inv[i+1] && w_inv[i] > w_inv[i+1]) return false;
    }
    
    return true;
}

bool is_bruhat_less_than(const std::vector<int>& x_data, const std::vector<int>& w_data) {
    if (x_data.size() != w_data.size()) return false;
    int n = x_data.size();
    if (n == 0) return true;

    std::vector<int> diff_row(n, 0);

    for (int i = 0; i < n; ++i) {
        int xi = x_data[i];
        int wi = w_data[i];
        
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