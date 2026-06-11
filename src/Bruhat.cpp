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
 * @brief Helper function to swap two int8_t values efficiently.
 * 
 * @param a Reference to the first value.
 * @param b Reference to the second value.
 */
inline void swap_val(ValueType& a, ValueType& b) {
    ValueType tmp = a;
    a = b;
    b = tmp;
}

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
bool make_ep_inplace(ValueType* x_data, ValueType* w_data, int n) {
    if (n == 0) return false;
    bool changed = true;
    bool changed_x = false;
    
    // Fixed-size stack buffers for calculating inverses (max size 32)
    ValueType x_inv[32];
    ValueType w_inv[32];
    
    // Compute inverses of the permutations
    for (int i = 0; i < n; ++i) {
        x_inv[x_data[i]] = i;
        w_inv[w_data[i]] = i;
    }
    
    while (changed) {
        changed = false;
        
        // Forward pass
        for (int i = 0; i < n - 1; ++i) {
            if (w_data[i] > w_data[i + 1] && x_data[i] < x_data[i + 1]) {
                int val1 = x_data[i];
                int val2 = x_data[i + 1];
                swap_val(x_data[i], x_data[i + 1]);
                swap_val(x_inv[val1], x_inv[val2]);
                changed = true;
            }
        }
        
        // Inverse pass
        for (int i = 0; i < n - 1; ++i) {
            if (w_inv[i] > w_inv[i + 1] && x_inv[i] < x_inv[i + 1]) {
                swap_val(x_data[x_inv[i]], x_data[x_inv[i + 1]]);
                swap_val(x_inv[i], x_inv[i + 1]);
                changed = true;
            }
        }
        
        if (changed) changed_x = true;
    }
    return changed_x;
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
    if (n == 0) return 0;

    // Frequency/difference array tracking value coverage (range 0..31)
    int dxw_row[32] = {0}; 
    int removed_count = 0;

    for (int i = 0; i < n; ++i) {
        if (x_data[i] < w_data[i]) {
            for (int j = x_data[i]; j < w_data[i]; ++j) {
                dxw_row[j]++;
            }
        } else {
            for (int j = w_data[i]; j < x_data[i]; ++j) {
                dxw_row[j]--;
            }
        }
        
        // Elements that are identical and have no overlapping intervals can be removed
        if (x_data[i] == w_data[i]) {
            if (dxw_row[x_data[i]] == 0) {
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

    // Map to keep track of active values
    int keep[32] = {0}; 

    for(int i = 0; i < n; ++i) {
        int v = x_data[i];
        if (v >= 0 && v < 32) keep[v] = 1; 
    }

    int value_map[32];
    std::fill(std::begin(value_map), std::end(value_map), REMOVED);

    int next = 0;
    // Build a re-indexing map for remaining valid values
    for (int val = 0; val < n; ++val) {
        if (keep[val]) value_map[val] = next++;
    }

    int write_idx = 0;
    // Compacting pass
    for (int read_idx = 0; read_idx < n; ++read_idx) {
        int xv = x_data[read_idx];
        if (xv != REMOVED) {
            x_data[write_idx] = static_cast<ValueType>(value_map[xv]);
            w_data[write_idx] = static_cast<ValueType>(value_map[w_data[read_idx]]);
            ++write_idx;
        }
    }
    n = write_idx; // Update physical size
}

} // namespace

// --- Public API Functions ---

void reduce(Permutation& x, Permutation& w) {
    auto* x_data = x.data_ptr_mut();
    auto* w_data = w.data_ptr_mut();
    int n = x.size();

    bool changed = true;
    bool changed_x = false;
    bool changed_w = false;
    
    while (changed) {
        if (make_ep_inplace(x_data, w_data, n)) changed_x = true;
        
        int removed_count = tilde_inplace(x_data, w_data, n);
        if (removed_count > 0) {
            flatten_inplace(x_data, w_data, n);
            // Updating size in wrapper Permutation objects
            x.set_size(n);
            w.set_size(n);
            changed = true;
            changed_x = true;
            changed_w = true;
        } else {
            changed = false;
        }
    }
    
    // Invalidate length caches if permutations were modified
    if (changed_x) x.reset_length_cache();
    if (changed_w) w.reset_length_cache();
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
    
    // Optimized implementation for n <= 32 using fixed arrays
    // Tracks the difference of the rank matrices: r_x(i, j) - r_w(i, j)
    int8_t diff_row[32] = {0};

    const auto* x_ptr = x.begin();
    const auto* w_ptr = w.begin();

    for (int i = 0; i < n; ++i) {
        int xi = x_ptr[i];
        int wi = w_ptr[i];
        
        // Branchless update loop optimized for basic processor instructions
        for (int j = 0; j < n; ++j) {
           // diff(i, j) = diff(i-1, j) + (j >= xi) - (j >= wi)
           // Boolean comparisons natively cast to 0 or 1
           diff_row[j] += (j >= xi) - (j >= wi);
           
           if (diff_row[j] < 0) return false;
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
    // Legacy implementation for vector inputs
    if (x_data.size() != w_data.size()) return false;
    int n = x_data.size();
    if (n == 0) return true;

    std::vector<int> prev_row(n, 0);
    std::vector<int> current_row(n, 0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            int val = (i == 0) ? 0 : prev_row[j];
            if (x_data[i] <= j && j < w_data[i]) val++;
            if (w_data[i] <= j && j < x_data[i]) val--;
            if (val < 0) return false;
            current_row[j] = val;
        }
        prev_row = current_row;
    }

    return true;
}

} // namespace Bruhat