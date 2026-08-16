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
constexpr int MAX_N_OPT = 16;

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
inline void compute_rx_matrix_flat(const Permutation& p, int8_t* buffer, int n) {
    if (n == 0) return;
    
    // Process the first row (stride MAX_N_OPT = 16)
    int p0 = p[0];
    for (int j = 0; j <= p0; ++j) buffer[j] = 1;
    for (int j = p0 + 1; j < MAX_N_OPT; ++j) buffer[j] = 0;
    
    // Process subsequent rows based on previous row results
    for (int i = 1; i < n; ++i) {
        const int8_t* prev_row = buffer + (i - 1) * MAX_N_OPT;
        int8_t* curr_row = buffer + i * MAX_N_OPT;
        int val = p[i];
        for (int j = 0; j <= val; ++j) {
            curr_row[j] = static_cast<int8_t>(prev_row[j] + 1);
        }
        for (int j = val + 1; j < MAX_N_OPT; ++j) {
            curr_row[j] = prev_row[j];
        }
    }
}

#include <immintrin.h>

alignas(16) static const int8_t MASK_LE_TABLE[16][16] = {
    {1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0},
    {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
};

/**
 * @brief Directly searches for extremal pairs z for (x, ws) with transposition s_index.
 */
template <typename Callback>
inline void search_extremal_pairs(const Permutation& x, const Permutation& ws, int s_index, Callback&& callback) {
    int n = x.size();
    if (n > MAX_N_OPT) {
        throw std::runtime_error("Permutation size exceeds optimized limit (16)");
    }

    const auto* ws_ptr = ws.begin();

    // Precompute rank matrices aligned to 16 bytes
    alignas(16) int8_t x_rx_buf[MAX_N_OPT * MAX_N_OPT] = {0};
    alignas(16) int8_t ws_rx_buf[MAX_N_OPT * MAX_N_OPT] = {0};
    compute_rx_matrix_flat(x, x_rx_buf, n);
    compute_rx_matrix_flat(ws, ws_rx_buf, n);

    // For ws_rx_buf, pad elements >= n with 127 so SIMD comparison rz <= ws_rx is always true for padding
    for (int row = 0; row < n; ++row) {
        int8_t* r = ws_rx_buf + row * MAX_N_OPT;
        for (int j = n; j < MAX_N_OPT; ++j) {
            r[j] = 127;
        }
    }

    // Precompute right and left descent bitmasks for ws
    unsigned int ws_rds_mask = 0;
    unsigned int ws_lds_mask = 0;
    
    int8_t ws_inv[MAX_N_OPT];
    for (int i = 0; i < n; ++i) {
        ws_inv[static_cast<int>(ws_ptr[i])] = static_cast<int8_t>(i);
    }
    
    for (int i = 0; i < n - 1; ++i) {
        if (ws_ptr[i] > ws_ptr[i + 1]) {
            ws_rds_mask |= (1u << i);
        }
        if (ws_inv[i] > ws_inv[i + 1]) {
            ws_lds_mask |= (1u << i);
        }
    }

    // Augmented descent mask: z must have descents where ws has descents, and also at s_index
    const unsigned int augmented_rds_mask = ws_rds_mask | (s_index >= 0 ? (1u << s_index) : 0);

    int8_t z_path[MAX_N_OPT];

    auto dfs = [&](auto self, int depth, unsigned int used_mask, int current_length, __m128i col) -> void {
        if (depth == n) {
            Permutation z(z_path, n);
            z.set_cached_length(current_length);
            callback(z);
            return;
        }

        const __m128i ws_rx_k = _mm_load_si128(reinterpret_cast<const __m128i*>(ws_rx_buf + depth * MAX_N_OPT));
        const __m128i x_rx_k = _mm_load_si128(reinterpret_cast<const __m128i*>(x_rx_buf + depth * MAX_N_OPT));

        int max_allowed_i = n - 1;
        if (depth > 0 && (augmented_rds_mask & (1u << (depth - 1)))) {
            int prev = static_cast<int>(z_path[depth - 1]);
            if (prev == 0) return; // 0 cannot descend to anything smaller
            max_allowed_i = prev - 1;
        }

        unsigned int max_mask = (1u << (max_allowed_i + 1)) - 1;
        unsigned int rejected_lds = ws_lds_mask & (~used_mask >> 1);
        unsigned int avail_mask = ~used_mask & max_mask & ~rejected_lds;

        alignas(16) int8_t col_arr[16];
        _mm_store_si128(reinterpret_cast<__m128i*>(col_arr), col);

        while (avail_mask) {
            int i = std::countr_zero(avail_mask);
            avail_mask &= avail_mask - 1; // clear lowest set bit
            
            unsigned int mask_bit = (1u << i);

            // Compute next rank row: rz = col + MASK_LE[i]
            __m128i mask_i = _mm_load_si128(reinterpret_cast<const __m128i*>(MASK_LE_TABLE[i]));
            __m128i rz = _mm_add_epi8(col, mask_i);

            // Bruhat check: rz <= ws_rx_k (invalid if rz > ws_rx_k)
            __m128i bad_ws = _mm_cmpgt_epi8(rz, ws_rx_k);
            // Bruhat check: x_rx_k <= rz (invalid if x_rx_k > rz)
            __m128i bad_x = _mm_cmpgt_epi8(x_rx_k, rz);

            __m128i bad = _mm_or_si128(bad_ws, bad_x);
            if (!_mm_testz_si128(bad, bad)) {
                continue;
            }

            int inv_added = col_arr[i + 1];
            z_path[depth] = static_cast<int8_t>(i);
            
            self(self, depth + 1, used_mask | mask_bit, current_length + inv_added, rz);
        }
    };

    dfs(dfs, 0, 0, 0, _mm_setzero_si128());
}

} // namespace

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

    // Exploit inverse symmetry: P_{x, w}(q) = P_{x^{-1}, w^{-1}}(q)
    Permutation x_inv = x.inverse();
    Permutation w_inv = w.inverse();
    if (std::make_pair(w.pack(), x.pack()) > std::make_pair(w_inv.pack(), x_inv.pack())) {
        x = x_inv;
        w = w_inv;
    }

    // Check the cache
    const auto key = std::make_pair(x.pack(), w.pack());
    if (auto it = cache.find(key); it != cache.end()) {
        return it->second;
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
    Polynomial result = get_kl_recursive(xs, ws, false);

    result.add_shifted(p_x_ws, 1);

    // Subtract summation terms
    apply_correction_term(result, x, ws, s_index);
    
    // Store the result in the cache
    cache.insert({key, result});
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
    
    int ws_len = ws.length();
    search_extremal_pairs(x, ws, s_index, [&](const Permutation& z) {
        int len_diff = ws_len - z.length();
        if (len_diff > 0 && (len_diff % 2 != 0)) {
            perform_correction(poly, x, ws, z);
        }
    });
}

void KLPolynomialCalculator::find_covering_relations_and_correct(
    Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index) {
    
    int n = x.size();
    
    for (int j = 0; j < n - 1; ++j) {
        for (int k = j + 1; k < n; ++k) {
            if (j != s_index && j != s_index + 1 && k != s_index && k != s_index + 1) {
                continue;
            }
            if (ws[j] > ws[k]) {
                // Early check: after swapping j and k, will z have a descent at s_index?
                int val_s = (j == s_index) ? ws[k] : ((k == s_index) ? ws[j] : ws[s_index]);
                int val_s1 = (j == s_index + 1) ? ws[k] : ((k == s_index + 1) ? ws[j] : ws[s_index + 1]);
                if (val_s <= val_s1) continue;

                bool covering = true;
                for (int l = j + 1; l < k; ++l) {
                    if (ws[l] > ws[k] && ws[l] < ws[j]) {
                        covering = false;
                        break;
                    }
                }
                
                if (covering) {
                    // Create permutation z by applying the covering relation swap
                    Permutation z = ws; 
                    std::swap(z.data_ptr_mut()[j], z.data_ptr_mut()[k]);
                    z.set_cached_length(ws.length() - 1); 
                    
                    if (Bruhat::is_bruhat_less_than(x, z) &&
                        !Bruhat::is_extremal_pair(z, ws)) {
                        
                        Polynomial p_x_z = get_kl_recursive(x, z, true);
                        poly.subtract_shifted(p_x_z, 1, 1);
                    }
                }
            }
        }
    }
}

void KLPolynomialCalculator::perform_correction(Polynomial& poly, const Permutation& x,
                                                const Permutation& ws, const Permutation& z) {
    int len_diff = ws.length() - z.length();
    if (len_diff <= 0 || (len_diff % 2 == 0)) return;

    long long mu = (len_diff == 1) ? 1 : get_mu_recursive(z, ws, true);
    if (mu == 0) return;

    Polynomial p_x_z = get_kl_recursive(x, z, true);
    int shift = (len_diff + 1) / 2;
    poly.subtract_shifted(p_x_z, shift, mu);
}

//##############################################################################
//##
//## Mu Coefficient Logic
//##
//##############################################################################

long long KLPolynomialCalculator::get_mu_coefficient(const Permutation& z, const Permutation& ws) {
    int len_diff = ws.length() - z.length();
    if (len_diff <= 0 || len_diff % 2 == 0) return 0;
    if (len_diff == 1) return 1;

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

    // Exploit inverse symmetry: mu(x, w) = mu(x^{-1}, w^{-1})
    Permutation x_inv = x.inverse();
    Permutation w_inv = w.inverse();
    if (std::make_pair(w.pack(), x.pack()) > std::make_pair(w_inv.pack(), x_inv.pack())) {
        x = x_inv;
        w = w_inv;
    }

    const auto key = std::make_pair(x.pack(), w.pack());
    if (auto it = mu_cache.find(key); it != mu_cache.end()) {
        return it->second;
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
    
    int ws_len = ws.length();
    int x_len = x.length();
    search_extremal_pairs(x, ws, s_index, [&](const Permutation& z) {
        int z_len = z.length();
        int len_diff_z_ws = ws_len - z_len;
        int len_diff_x_z = z_len - x_len;
        if (len_diff_z_ws > 1 && (len_diff_z_ws % 2 != 0) &&
            len_diff_x_z > 0 && (len_diff_x_z % 2 != 0)) {
            mu_perform_correction(current_mu, x, ws, z);
        }
    });
}

void KLPolynomialCalculator::mu_find_covering_relations_and_correct(
    long long& current_mu, const Permutation& x, const Permutation& ws, int s_index) {
    
    int n = x.size();
    
    for (int j = 0; j < n - 1; ++j) {
        for (int k = j + 1; k < n; ++k) {
            if (j != s_index && j != s_index + 1 && k != s_index && k != s_index + 1) {
                continue;
            }
            if (ws[j] > ws[k]) {
                // Early check: after swapping j and k, will z have a descent at s_index?
                int val_s = (j == s_index) ? ws[k] : ((k == s_index) ? ws[j] : ws[s_index]);
                int val_s1 = (j == s_index + 1) ? ws[k] : ((k == s_index + 1) ? ws[j] : ws[s_index + 1]);
                if (val_s <= val_s1) continue;

                bool covering = true;
                for (int l = j + 1; l < k; ++l) {
                    if (ws[l] > ws[k] && ws[l] < ws[j]) {
                        covering = false;
                        break;
                    }
                }
                if (covering) {
                    // Fast element swap setup
                    Permutation z = ws;
                    std::swap(z.data_ptr_mut()[j], z.data_ptr_mut()[k]);
                    z.set_cached_length(ws.length() - 1);
                    
                    if (Bruhat::is_bruhat_less_than(x, z)) {
                        int len_diff_x_z = z.length() - x.length();
                        if (len_diff_x_z > 0 && (len_diff_x_z % 2 != 0)) {
                            long long mu_x_z = (len_diff_x_z == 1) ? 1 : get_mu_recursive(x, z, true);
                            current_mu -= mu_x_z;
                            if (current_mu == 0) return;
                        }
                    }
                }
            }
        }
    }
}

void KLPolynomialCalculator::mu_perform_correction(
    long long& current_mu, const Permutation& x, const Permutation& ws, const Permutation& z) {
    int len_diff_z_ws = ws.length() - z.length();
    if (len_diff_z_ws <= 0 || (len_diff_z_ws % 2 == 0)) return;

    int len_diff_x_z = z.length() - x.length();
    if (len_diff_x_z <= 0 || (len_diff_x_z % 2 == 0)) return;

    long long mu_z_ws = (len_diff_z_ws == 1) ? 1 : get_mu_recursive(z, ws, true);
    if (mu_z_ws == 0) return;

    long long mu_x_z = (len_diff_x_z == 1) ? 1 : get_mu_recursive(x, z, true);
    current_mu -= mu_z_ws * mu_x_z;
}

