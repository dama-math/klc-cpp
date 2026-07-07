#include "KLCalculator.hpp"
#include "Bruhat.hpp"
#include "Lascoux.hpp"
#include <algorithm> 
#include <numeric>   
#include <iostream>
#include <array>
#include <functional>
#include <bit>

namespace {

// Constant defining the maximum permutation size supported by optimized fixed-size routines
constexpr int MAX_N_OPT = 32;

/**
 * @brief Helper function to compute the rank matrix efficiently from Permutation iterators.
 * 
 * Computes the rank matrix for a given permutation and stores the result in a flat,
 * row-major array buffer to optimize cache locality.
 * 
 * @param p The permutation to process.
 * @param buffer The destination buffer (must be of size n * n).
 * @param n The size of the permutation.
 */
void compute_rx_matrix_flat(const Permutation& p, int* buffer, int n) {
    if (n == 0) return;
    
    // Process the first row
    int p0 = p[0];
    for (int j = 0; j < n; ++j) {
        buffer[j] = (p0 >= j) ? 1 : 0;
    }
    
    // Process subsequent rows based on previous row results
    for (int i = 1; i < n; ++i) {
        const int* prev_row = buffer + (i - 1) * n;
        int* curr_row = buffer + i * n;
        int val = p[i];
        for (int j = 0; j < n; ++j) {
            curr_row[j] = prev_row[j] + ((val >= j) ? 1 : 0);
        }
    }
}

/**
 * @brief Encapsulates the DFS state for finding extremal pairs during correction term derivation.
 * 
 * This class abstracts the highly optimized depth-first search (DFS) procedure 
 * utilized to find valid permutations `z` that form extremal pairs with `ws`.
 */
class ExtremalPairSearcher {
public:
    /**
     * @brief Constructs an ExtremalPairSearcher and precomputes required matrices.
     * 
     * @param x Inner permutation constraint.
     * @param ws The modified outer permutation.
     * @param s_index The index of the simple transposition applied.
     */
    ExtremalPairSearcher(const Permutation& x, const Permutation& ws, int s_index)
        : n(x.size()), s_index(s_index) {
        
        if (n > MAX_N_OPT) {
            throw std::runtime_error("Permutation size exceeds optimized limit (32)");
        }

        // Precompute rank matrices
        compute_rx_matrix_flat(x, x_rx_buf.data(), n);
        compute_rx_matrix_flat(ws, ws_rx_buf.data(), n);

        // Precompute right and left descent information for ws
        ws_rds.fill(false);
        ws_lds.fill(false);
        ws_lds_mask = 0;
        
        Permutation ws_inv = ws.inverse();
        
        for (int i = 0; i < n - 1; ++i) {
            if (ws[i] > ws[i + 1]) ws_rds[i] = true;
            if (ws_inv[i] > ws_inv[i + 1]) {
                ws_lds[i] = true;
                ws_lds_mask |= (1u << i);
            }
        }

        // Initialize row pointers for fast direct memory access
        for (int r = 0; r < n; ++r) {
            x_rx_ptrs[r] = x_rx_buf.data() + r * n;
            ws_rx_ptrs[r] = ws_rx_buf.data() + r * n;
        }
        
        // Copy permutation data to fast local arrays
        for (int i = 0; i < n; ++i) {
            x_data[i] = x[i];
            ws_data[i] = ws[i];
        }
    }

    /**
     * @brief Initiates the DFS search and applies the callback to found permutations.
     * 
     * @tparam Callback Type of the callback function.
     * @param callback Function to execute when a valid permutation `z` is found.
     */
    template <typename Callback>
    void search(Callback callback) {
        std::array<int8_t, MAX_N_OPT> curr_col{};
        curr_col.fill(0); 
        
        std::array<int8_t, MAX_N_OPT> z_path;
        
        dfs(0, 0, curr_col, z_path, callback);
    }

private:
    int n;
    int s_index;

    // Fixed-size arrays to avoid dynamic heap allocation during intense recursion
    std::array<int, MAX_N_OPT * MAX_N_OPT> x_rx_buf;
    std::array<int, MAX_N_OPT * MAX_N_OPT> ws_rx_buf;
    std::array<bool, MAX_N_OPT> ws_rds;
    std::array<bool, MAX_N_OPT> ws_lds;
    
    // Fast access pointers/arrays
    const int* x_rx_ptrs[MAX_N_OPT]{};
    const int* ws_rx_ptrs[MAX_N_OPT]{};
    int x_data[MAX_N_OPT]{};
    int ws_data[MAX_N_OPT]{};
    unsigned int ws_lds_mask = 0;

    /**
     * @brief The core recursive DFS function.
     */
    template <typename Callback>
    void dfs(int depth, unsigned int used_mask, 
             std::array<int8_t, MAX_N_OPT>& prev_col, 
             std::array<int8_t, MAX_N_OPT>& z_path,
             Callback& callback) {
        
        if (depth == n) {
            // Reconstruct the permutation from the finalized z_path
            Permutation z(z_path.data(), n);
            callback(z);
            return;
        }

        const int ws_k = ws_data[depth];
        const int x_k = x_data[depth];
        const int* ws_rx_k = ws_rx_ptrs[depth];
        const int* x_rx_k = x_rx_ptrs[depth];

        int max_allowed_i = n - 1;
        if (depth > 0 && ws_rds[depth - 1]) {
            max_allowed_i = std::min(max_allowed_i, static_cast<int>(z_path[depth - 1]));
        }
        if (s_index != -1 && depth == s_index + 1) {
            max_allowed_i = std::min(max_allowed_i, static_cast<int>(z_path[s_index]));
        }

        unsigned int max_mask = (1u << (max_allowed_i + 1)) - 1;
        unsigned int rejected_lds = ws_lds_mask & (~used_mask >> 1);
        unsigned int avail_mask = ~used_mask & max_mask & ~rejected_lds;

        while (avail_mask) {
            int i = std::countr_zero(avail_mask);
            avail_mask &= avail_mask - 1; // clear lowest set bit
            
            unsigned int mask_bit = (1u << i);

            bool bruhat_ok = true;
            
            // Check Bruhat constraints against ws
            for (int j = ws_k + 1; j <= i; ++j) {
                if (ws_rx_k[j] <= prev_col[j]) { bruhat_ok = false; break; }
            }
            if (!bruhat_ok) continue;

            // Check Bruhat constraints against x
            for (int j = i + 1; j <= x_k; ++j) {
                if (prev_col[j] < x_rx_k[j]) { bruhat_ok = false; break; }
            }
            if (!bruhat_ok) continue;

            // Update column counts and recurse
            for (int k = 0; k <= i; ++k) prev_col[k]++;
            z_path[depth] = static_cast<int8_t>(i);
            
            dfs(depth + 1, used_mask | mask_bit, prev_col, z_path, callback);
            
            // Undo column counts
            for (int k = 0; k <= i; ++k) prev_col[k]--;
        }
    }
};

} // namespace

//##############################################################################
//##
//## Public API Implementation
//##
//##############################################################################

KLPolynomialCalculator::KLPolynomialCalculator() = default;

Polynomial KLPolynomialCalculator::calculate_kl(const Permutation& x, const Permutation& w, bool use_lascoux) {
    if (x.size() != w.size()) {
        throw std::invalid_argument("Permutations x and w must have the same size.");
    }
    
    // Polynomial is zero if x is not less than or equal to w in the Bruhat order
    if (!Bruhat::is_bruhat_less_than(x, w)) {
        return Polynomial();
    }

    // Polynomial is 1 if the length difference is 0, 1, or 2
    if (w.length() - x.length() <= 2) {
        return Polynomial({1});
    }

    // Use Lascoux's optimization if enabled, applicable, and w is 3412-avoiding
    if (use_lascoux && x.size() >= 8 && !Lascoux::has_3412_pattern(w)) {
        return Lascoux::calculate_kl(x, w);
    } else {
        return get_kl_recursive(x, w, true);
    }
}

void KLPolynomialCalculator::clear_cache() {
    cache.clear();
    mu_cache.clear();
    poly_pool.clear();
}

void KLPolynomialCalculator::precompute_kl_polynomials(int max_n) {
    for (int n = 1; n <= max_n; ++n) {
        std::vector<int> p(n);
        std::iota(p.begin(), p.end(), 0);

        std::vector<Permutation> permutations;
        do {
            permutations.emplace_back(p);
        } while (std::next_permutation(p.begin(), p.end()));

        for (const auto& w : permutations) {
            for (const auto& x : permutations) {
                if (Bruhat::is_bruhat_less_than(x, w)) {
                    calculate_kl(x, w);
                }
            }
        }
    }
}

//##############################################################################
//##
//## Main Recursive Logic
//##
//##############################################################################

Polynomial KLPolynomialCalculator::get_kl_recursive(Permutation x, Permutation w, bool already_bruhat_checked) {
    if (!already_bruhat_checked && !Bruhat::is_bruhat_less_than(x, w)) {
        return Polynomial();
    }

    if (w.length() - x.length() <= 2) {
        return Polynomial({1});
    }

    // Attempt to reduce the pair to a simpler equivalent pair
    Bruhat::reduce(x, w);

    if (w.length() - x.length() <= 2) {
        return Polynomial({1});
    }

    // Check the cache
    const auto key = std::make_pair(x.pack(), w.pack());
    if (const auto* cached_id = cache.find(key)) {
        return poly_pool.get_poly(*cached_id);
    }

    int s_index  = -1;
    int pos_zero = -1;
    int pos_max  = -1;

    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0)             pos_zero = i;
        if (w[i] == w.size() - 1)  pos_max  = i;
    }
    
    if (pos_zero > 0) {
        s_index = (pos_max > w.size() - 1 - pos_zero) ? pos_max : (pos_zero - 1);
    }

    if (s_index == -1) {
        return Polynomial({1}); // Should not happen for length diff > 2, but safe fallback
    }

    Permutation ws = w.multiply_simple_transposition(s_index);
    Permutation xs = x.multiply_simple_transposition(s_index);

    // Core recurrence relation calculation
    Polynomial p_x_ws = get_kl_recursive(x, ws, false);
    Polynomial p_xs_ws = get_kl_recursive(xs, ws, false);

    Polynomial result = p_xs_ws;
    result += p_x_ws.shift(1); 

    // Subtract summation terms
    apply_correction_term(result, x, ws, s_index);
    
    // Store the result in the cache and polynomial pool
    uint32_t result_id = poly_pool.get_id(result);
    cache[key] = result_id;
    mu_cache[key] = result.get_coeff((w.length() - x.length() - 1) / 2);
    return result;
}

//##############################################################################
//##
//## Correction Term Logic
//##
//##############################################################################

void KLPolynomialCalculator::apply_correction_term(Polynomial& poly, const Permutation& x,
                                                   const Permutation& ws, int s_index) {
    int n = x.size();
    if (n == 0) return;

    find_extremal_pairs_and_correct(poly, x, ws, s_index);
    find_covering_relations_and_correct(poly, x, ws, s_index);
}

void KLPolynomialCalculator::find_extremal_pairs_and_correct(
    Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index) {
    
    ExtremalPairSearcher searcher(x, ws, s_index);
    searcher.search([&](Permutation z) {
        perform_correction(poly, x, ws, z);
    });
}

void KLPolynomialCalculator::find_covering_relations_and_correct(
    Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index) {
    
    int n = x.size();
    
    // Search for covering relations by checking possible swaps
    for (int j = 0; j < n - 1; ++j) {
        for (int k = j + 1; k < n; ++k) {
            if (ws[j] > ws[k]) {
                int inversions_between = 0;
                for (int l = j + 1; l < k; ++l) {
                    if (ws[l] > ws[k] && ws[l] < ws[j]) {
                        inversions_between++;
                    }
                }
                
                if (inversions_between == 0) {
                    // Create permutation z by applying the covering relation swap
                    Permutation z = ws; 
                    std::swap(z.data_ptr_mut()[j], z.data_ptr_mut()[k]);
                    z.reset_length_cache(); 
                    
                    if (z[s_index] > z[s_index + 1] &&
                        Bruhat::is_bruhat_less_than(x, z) &&
                        !Bruhat::is_extremal_pair(z, ws)) {
                        
                        perform_correction(poly, x, ws, z);
                    }
                }
            }
        }
    }
}

void KLPolynomialCalculator::perform_correction(Polynomial& poly, const Permutation& x,
                                                const Permutation& ws, const Permutation& z) {
    long long mu = get_mu_recursive(z, ws, true);
    if (mu == 0) return;

    Polynomial p_x_z = get_kl_recursive(x, z, true);
    int shift = (ws.length() - z.length() + 1) / 2;
    Polynomial correction = p_x_z.shift(shift).multiply_by_scalar(mu);
    poly -= correction;
}

//##############################################################################
//##
//## Mu Coefficient Logic
//##
//##############################################################################

long long KLPolynomialCalculator::get_mu_coefficient(const Permutation& z, const Permutation& ws) {
    int len_diff = ws.length() - z.length();
    if (len_diff <= 0 || len_diff % 2 == 0) return 0;

    Polynomial p_z_ws = get_kl_recursive(z, ws, true);
    int deg = p_z_ws.degree();

    if (deg >= 0 && (2 * deg == len_diff - 1)) {
        return p_z_ws.get_coeff(deg);
    }
    return 0;
}

long long KLPolynomialCalculator::get_mu_recursive(Permutation x, Permutation w, bool already_bruhat_checked) {
    if (!already_bruhat_checked && !Bruhat::is_bruhat_less_than(x, w)) {
        return 0;
    }

    int len_diff = w.length() - x.length();
    if (len_diff <= 0 || len_diff % 2 == 0) return 0;
    if (len_diff == 1) return 1;

    // Apply pair reduction and check if originally an extremal pair
    if (!Bruhat::reduce(x, w)) return 0;

    int reduced_len_diff = w.length() - x.length();
    if (reduced_len_diff <= 0 || reduced_len_diff % 2 == 0) return 0;
    if (reduced_len_diff == 1) return 1;

    if (reduced_len_diff < len_diff) {
        return 0; 
    }

    const auto key = std::make_pair(x.pack(), w.pack());
    if (const auto* cached_mu = mu_cache.find(key)) {
        return *cached_mu;
    }

    int s_index  = -1;
    int pos_zero = -1;
    int pos_max  = -1;

    for (int i = 0; i < w.size(); ++i) {
        if (w[i] == 0)             pos_zero = i;
        if (w[i] == w.size() - 1)  pos_max  = i;
    }

    if (pos_zero > 0) {
        s_index = (pos_max > w.size() - 1 - pos_zero) ? pos_max : (pos_zero - 1);
    }

    if (s_index == -1) {
        return 1;
    }

    Permutation ws = w.multiply_simple_transposition(s_index);
    Permutation xs = x.multiply_simple_transposition(s_index);

    Polynomial p_x_ws = get_kl_recursive(x, ws, false);
    long long mu_xs_ws = get_mu_recursive(xs, ws, false);

    long long result = p_x_ws.get_coeff((len_diff - 3) / 2);
    result += mu_xs_ws;

    if (result == 0) {
        mu_cache[key] = 0;
        return 0;
    }

    // Apply corrections specific to mu derivation
    mu_apply_correction_term(result, x, ws, s_index);

    mu_cache[key] = result;
    return result;
}

void KLPolynomialCalculator::mu_apply_correction_term(long long& current_mu, const Permutation& x,
                                                   const Permutation& ws, int s_index) {
    int n = x.size();
    if (n == 0) return;

    mu_find_covering_relations_and_correct(current_mu, x, ws, s_index);

    if (current_mu == 0) return; 

    mu_find_extremal_pairs_and_correct(current_mu, x, ws, s_index);
}

void KLPolynomialCalculator::mu_find_extremal_pairs_and_correct(
    long long& current_mu, const Permutation& x, const Permutation& ws, int s_index) {
    
    ExtremalPairSearcher searcher(x, ws, s_index);
    searcher.search([&](Permutation z) {
        if (z.length() < ws.length() - 1) { 
             mu_perform_correction(current_mu, x, ws, z);
        }
    });
}

void KLPolynomialCalculator::mu_find_covering_relations_and_correct(
    long long& current_mu, const Permutation& x, const Permutation& ws, int s_index) {
    
    int n = x.size();
    
    for (int j = 0; j < n - 1; ++j) {
        for (int k = j + 1; k < n; ++k) {
            if (ws[j] > ws[k]) {
                int inversions_between = 0;
                for (int l = j + 1; l < k; ++l) {
                    if (ws[l] > ws[k] && ws[l] < ws[j]) {
                        inversions_between++;
                    }
                }
                if (inversions_between == 0) {
                    // Fast element swap setup
                    Permutation z = ws;
                    std::swap(z.data_ptr_mut()[j], z.data_ptr_mut()[k]);
                    z.reset_length_cache();
                    
                    if (z[s_index] > z[s_index + 1] && Bruhat::is_bruhat_less_than(x, z)) {
                        mu_perform_correction(current_mu, x, ws, z);
                        if (current_mu == 0) return;
                    }
                }
            }
        }
    }
}

void KLPolynomialCalculator::mu_perform_correction(
    long long& current_mu, const Permutation& x, const Permutation& ws, const Permutation& z) {
    long long mu_z_ws = get_mu_recursive(z, ws, true);
    if (mu_z_ws == 0) return;

    long long mu_x_z = get_mu_recursive(x, z, true);
    current_mu -= mu_z_ws * mu_x_z;
}

std::vector<std::vector<int>> KLPolynomialCalculator::compute_rx_matrix(const std::vector<int>& p_data) {
    int n = p_data.size();
    std::vector<int> buf(n * n);
    
    // Utilize the flat optimized routine under the hood
    Permutation p(p_data);
    compute_rx_matrix_flat(p, buf.data(), n);
    
    std::vector<std::vector<int>> res(n, std::vector<int>(n));
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            res[i][j] = buf[i * n + j];
        }
    }
    return res;
}
