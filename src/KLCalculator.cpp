#include "KLCalculator.hpp"
#include "Bruhat.hpp"
#include "Lascoux.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <immintrin.h>
#include <stdexcept>

namespace {

constexpr int MAX_N_OPT = 16;

/**
 * @brief Computes the 2D rank matrix rank(p[1..i, 1..j]) into a row-major flat buffer.
 * Stride between consecutive rows is MAX_N_OPT (16 bytes) for direct SIMD loads.
 */
// 16-byte aligned indicator vectors: row i has 1s in columns 0..i and 0s elsewhere
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
 * @brief Computes the 2D rank matrix rank(p[1..i, 1..j]) into a row-major flat buffer.
 * Stride between consecutive rows is MAX_N_OPT (16 bytes) for direct SIMD loads.
 */
inline void compute_rx_matrix_flat(const Permutation& p, int8_t* buffer, int n) {
    if (n == 0) return;
    
    __m128i row = _mm_load_si128(reinterpret_cast<const __m128i*>(MASK_LE_TABLE[p[0]]));
    _mm_store_si128(reinterpret_cast<__m128i*>(buffer), row);
    
    for (int i = 1; i < n; ++i) {
        __m128i mask = _mm_load_si128(reinterpret_cast<const __m128i*>(MASK_LE_TABLE[p[i]]));
        row = _mm_add_epi8(row, mask);
        _mm_store_si128(reinterpret_cast<__m128i*>(buffer + i * MAX_N_OPT), row);
    }
}

/**
 * @brief SIMD depth-first search for permutations z in [x, ws] that are extremal w.r.t. ws
 * and have a descent at s_index.
 */
template <typename Callback>
inline void search_extremal_pairs(const Permutation& x, const Permutation& ws, int s_index, Callback&& callback) {
    int n = x.size();
    if (n > MAX_N_OPT) {
        throw std::runtime_error("Permutation size exceeds optimized limit (16)");
    }

    const auto* ws_ptr = ws.begin();

    alignas(16) int8_t x_rx_buf[MAX_N_OPT * MAX_N_OPT];
    alignas(16) int8_t ws_rx_buf[MAX_N_OPT * MAX_N_OPT];
    compute_rx_matrix_flat(x, x_rx_buf, n);
    compute_rx_matrix_flat(ws, ws_rx_buf, n);

    // Pad elements beyond n with 127 so SIMD comparison rz <= ws_rx is always true for padding
    for (int row = 0; row < n; ++row) {
        int8_t* r = ws_rx_buf + row * MAX_N_OPT;
        for (int j = n; j < MAX_N_OPT; ++j) {
            r[j] = 127;
        }
    }

    // Precompute right and left descent bitmasks for ws using SIMD
    __m128i ws_v = _mm_loadu_si128(reinterpret_cast<const __m128i*>(ws_ptr));
    __m128i ws_v_shift = _mm_srli_si128(ws_v, 1);
    unsigned int ws_rds_mask = _mm_movemask_epi8(_mm_cmpgt_epi8(ws_v, ws_v_shift)) & ((1u << (n - 1)) - 1);

    alignas(16) int8_t ws_inv[MAX_N_OPT];
    for (int i = 0; i < n; ++i) {
        ws_inv[static_cast<int>(ws_ptr[i])] = static_cast<int8_t>(i);
    }
    __m128i inv_v = _mm_load_si128(reinterpret_cast<const __m128i*>(ws_inv));
    __m128i inv_v_shift = _mm_srli_si128(inv_v, 1);
    unsigned int ws_lds_mask = _mm_movemask_epi8(_mm_cmpgt_epi8(inv_v, inv_v_shift)) & ((1u << (n - 1)) - 1);

    // Augmented descent mask: z must have descents where ws has descents, plus at s_index
    const unsigned int augmented_rds_mask = ws_rds_mask | (s_index >= 0 ? (1u << s_index) : 0);
    int8_t z_path[MAX_N_OPT];

    if (n == 0) return;
    if (n == 1) {
        Permutation z = x;
        callback(z);
        return;
    }

    auto dfs = [&](auto self, int depth, unsigned int used_mask, int current_length, __m128i col) -> bool {
        const __m128i ws_rx_k = _mm_load_si128(reinterpret_cast<const __m128i*>(ws_rx_buf + depth * MAX_N_OPT));
        const __m128i x_rx_k = _mm_load_si128(reinterpret_cast<const __m128i*>(x_rx_buf + depth * MAX_N_OPT));

        int max_allowed_i = n - 1;
        if (depth > 0 && (augmented_rds_mask & (1u << (depth - 1)))) {
            int prev = static_cast<int>(z_path[depth - 1]);
            if (prev == 0) return true;
            max_allowed_i = prev - 1;
        }

        unsigned int max_mask = (1u << (max_allowed_i + 1)) - 1;
        unsigned int rejected_lds = ws_lds_mask & (~used_mask >> 1);
        unsigned int avail_mask = ~used_mask & max_mask & ~rejected_lds;

        if (depth == n - 1) {
            if (avail_mask) {
                int i = std::countr_zero(avail_mask);
                __m128i mask_i = _mm_load_si128(reinterpret_cast<const __m128i*>(MASK_LE_TABLE[i]));
                __m128i rz = _mm_add_epi8(col, mask_i);

                __m128i bad_ws = _mm_cmpgt_epi8(rz, ws_rx_k);
                __m128i bad_x = _mm_cmpgt_epi8(x_rx_k, rz);
                __m128i bad = _mm_or_si128(bad_ws, bad_x);
                if (_mm_testz_si128(bad, bad)) {
                    int inv_added = std::popcount(used_mask >> (i + 1));
                    z_path[n - 1] = static_cast<int8_t>(i);
                    Permutation z(z_path, n);
                    z.set_cached_length(current_length + inv_added);
                    return callback(z);
                }
            }
            return true;
        }

        while (avail_mask) {
            int i = std::countr_zero(avail_mask);
            avail_mask &= avail_mask - 1;
            
            unsigned int mask_bit = (1u << i);

            // Compute rank row vector for candidate: rz = col + indicator(<= i)
            __m128i mask_i = _mm_load_si128(reinterpret_cast<const __m128i*>(MASK_LE_TABLE[i]));
            __m128i rz = _mm_add_epi8(col, mask_i);

            // Verify Bruhat rank condition: x_rx_k <= rz <= ws_rx_k
            __m128i bad_ws = _mm_cmpgt_epi8(rz, ws_rx_k);
            __m128i bad_x = _mm_cmpgt_epi8(x_rx_k, rz);
            __m128i bad = _mm_or_si128(bad_ws, bad_x);
            if (!_mm_testz_si128(bad, bad)) {
                continue;
            }

            int inv_added = std::popcount(used_mask >> (i + 1));
            z_path[depth] = static_cast<int8_t>(i);
            
            if (!self(self, depth + 1, used_mask | mask_bit, current_length + inv_added, rz)) {
                return false;
            }
        }
        return true;
    };

    dfs(dfs, 0, 0, 0, _mm_setzero_si128());
}

/**
 * @brief Selects the optimal simple reflection s_index = (i, i+1) where w has a descent.
 * Warrington's heuristic picks s near 0 or n-1 to maximize subproblem reuse.
 */
inline int select_descent_index(const Permutation& w) {
    int pos_zero = -1;
    int pos_max = -1;
    int n = w.size();

    for (int i = 0; i < n; ++i) {
        if (w[i] == 0)     pos_zero = i;
        if (w[i] == n - 1) pos_max  = i;
    }

    if (pos_zero > 0) {
        return (pos_max > n - 1 - pos_zero) ? pos_max : (pos_zero - 1);
    }
    return -1;
}

/**
 * @brief Enumerates all covering relations z <| ws where z has a descent at s_index.
 * Invokes callback(z). If callback returns false, stops traversal early.
 */
template <typename Callback>
inline void for_each_covering_descent(const Permutation& ws, int s_index, Callback&& callback) {
    int n = ws.size();
    
    for (int j = 0; j < n - 1; ++j) {
        for (int k = j + 1; k < n; ++k) {
            // Must touch the transposition interval {s_index, s_index + 1}
            if (j != s_index && j != s_index + 1 && k != s_index && k != s_index + 1) {
                continue;
            }
            if (ws[j] > ws[k]) {
                // Check if swapping (j, k) produces a descent at s_index
                int val_s  = (j == s_index) ? ws[k] : ((k == s_index) ? ws[j] : ws[s_index]);
                int val_s1 = (j == s_index + 1) ? ws[k] : ((k == s_index + 1) ? ws[j] : ws[s_index + 1]);
                if (val_s <= val_s1) continue;

                // Check covering condition: no intermediate element between ws[k] and ws[j]
                bool covering = true;
                for (int l = j + 1; l < k; ++l) {
                    if (ws[l] > ws[k] && ws[l] < ws[j]) {
                        covering = false;
                        break;
                    }
                }

                if (covering) {
                    Permutation z = ws;
                    std::swap(z.data_ptr_mut()[j], z.data_ptr_mut()[k]);
                    z.set_cached_length(ws.length() - 1);

                    if (!callback(z)) return;
                }
            }
        }
    }
}

} // anonymous namespace

KLPolynomialCalculator::KLPolynomialCalculator() = default;

Polynomial KLPolynomialCalculator::calculate_kl(const Permutation& x, const Permutation& w, bool use_lascoux) {
    if (x.size() != w.size()) {
        throw std::invalid_argument("Permutations x and w must have the same size.");
    }
    
    if (!Bruhat::is_bruhat_less_than(x, w)) {
        return Polynomial();
    }

    if (w.length() - x.length() <= 2) {
        return Polynomial({1});
    }

    // Use Lascoux's fast method for vexillary (3412-avoiding) permutations when n >= 8
    if (use_lascoux && x.size() >= 8 && !Lascoux::has_3412_pattern(w)) {
        return Lascoux::calculate_kl(x, w);
    }

    return get_kl_recursive(x, w, true);
}

void KLPolynomialCalculator::clear_cache() {
    cache.clear();
    mu_cache.clear();
}

//##############################################################################
//## Core Recursive Solver for P_{x, w}(q)
//##############################################################################

Polynomial KLPolynomialCalculator::get_kl_recursive(Permutation x, Permutation w, bool already_bruhat_checked) {
    if (!already_bruhat_checked && !Bruhat::is_bruhat_less_than(x, w)) {
        return Polynomial();
    }

    if (w.length() - x.length() <= 2) {
        return Polynomial({1});
    }

    // Apply extremal pair and tilde reductions
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

    const auto key = std::make_pair(x.pack(), w.pack());
    if (auto* cached = cache.find(key)) {
        return *cached;
    }

    int s_index = select_descent_index(w);
    if (s_index == -1) {
        return Polynomial({1});
    }

    Permutation ws = w.multiply_simple_transposition(s_index);
    Permutation xs = x.multiply_simple_transposition(s_index);

    // Recurrence: P_{x, w} = q * P_{x, ws} + P_{xs, ws} - corrections
    Polynomial p_x_ws = get_kl_recursive(x, ws, false);
    Polynomial result = get_kl_recursive(xs, ws, false);
    result.add_shifted(p_x_ws, 1);

    apply_correction_term(result, x, ws, s_index);
    
    cache.insert_or_assign(key, result);
    mu_cache.insert_or_assign(key, result.get_coeff((w.length() - x.length() - 1) / 2));
    return result;
}

void KLPolynomialCalculator::apply_correction_term(Polynomial& poly, const Permutation& x,
                                                   const Permutation& ws, int s_index) {
    if (x.size() == 0) return;
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
        return true;
    });
}

void KLPolynomialCalculator::find_covering_relations_and_correct(
    Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index) {
    
    for_each_covering_descent(ws, s_index, [&](const Permutation& z) {
        if (Bruhat::is_bruhat_less_than(x, z) && !Bruhat::is_extremal_pair(z, ws)) {
            Polynomial p_x_z = get_kl_recursive(x, z, true);
            poly.subtract_shifted(p_x_z, 1, 1);
        }
        return true;
    });
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
//## Core Recursive Solver for mu(x, w)
//##############################################################################

long long KLPolynomialCalculator::get_mu_recursive(Permutation x, Permutation w, bool already_bruhat_checked) {
    if (!already_bruhat_checked && !Bruhat::is_bruhat_less_than(x, w)) {
        return 0;
    }

    int len_diff = w.length() - x.length();
    if (len_diff <= 0 || len_diff % 2 == 0) return 0;
    if (len_diff == 1) return 1;

    // Apply pair reduction: if not originally an extremal pair, mu is 0
    if (!Bruhat::reduce(x, w)) return 0;

    int reduced_len_diff = w.length() - x.length();
    if (reduced_len_diff <= 0 || reduced_len_diff % 2 == 0) return 0;
    if (reduced_len_diff == 1) return 1;
    if (reduced_len_diff < len_diff) return 0;

    // Exploit inverse symmetry: mu(x, w) = mu(x^{-1}, w^{-1})
    Permutation x_inv = x.inverse();
    Permutation w_inv = w.inverse();
    if (std::make_pair(w.pack(), x.pack()) > std::make_pair(w_inv.pack(), x_inv.pack())) {
        x = x_inv;
        w = w_inv;
    }

    const auto key = std::make_pair(x.pack(), w.pack());
    if (auto* cached = mu_cache.find(key)) {
        return *cached;
    }

    int s_index = select_descent_index(w);
    if (s_index == -1) {
        return 1;
    }

    Permutation ws = w.multiply_simple_transposition(s_index);
    Permutation xs = x.multiply_simple_transposition(s_index);

    Polynomial p_x_ws = get_kl_recursive(x, ws, false);
    long long mu_xs_ws = get_mu_recursive(xs, ws, false);

    long long result = p_x_ws.get_coeff((len_diff - 3) / 2) + mu_xs_ws;
    if (result == 0) {
        mu_cache.insert_or_assign(key, 0);
        return 0;
    }

    mu_apply_correction_term(result, x, ws, s_index);
    mu_cache.insert_or_assign(key, result);
    return result;
}

void KLPolynomialCalculator::mu_apply_correction_term(long long& current_mu, const Permutation& x,
                                                      const Permutation& ws, int s_index) {
    if (x.size() == 0) return;

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
            if (current_mu == 0) return false;
        }
        return true;
    });
}

void KLPolynomialCalculator::mu_find_covering_relations_and_correct(
    long long& current_mu, const Permutation& x, const Permutation& ws, int s_index) {
    
    for_each_covering_descent(ws, s_index, [&](const Permutation& z) {
        if (Bruhat::is_bruhat_less_than(x, z)) {
            int len_diff_x_z = z.length() - x.length();
            if (len_diff_x_z > 0 && (len_diff_x_z % 2 != 0)) {
                long long mu_x_z = (len_diff_x_z == 1) ? 1 : get_mu_recursive(x, z, true);
                current_mu -= mu_x_z;
                if (current_mu == 0) return false;
            }
        }
        return true;
    });
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
