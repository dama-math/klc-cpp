#include "KLCalculator.hpp"
#include "Bruhat.hpp"
#include "Lascoux.hpp"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include "detail/IntervalSearch.hpp"
#include "detail/CheckedArithmetic.hpp"
#include <stdexcept>

namespace {

using klc::detail::search_extremal_pairs;

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
                    Permutation z = klc::detail::PermutationAccess::covering_swap(ws, j, k);

                    if (!callback(z)) return;
                }
            }
        }
    }
}

// Both recurrences use the same inverse symmetry and key ordering.
inline auto canonicalize_pair(Permutation& x, Permutation& w) {
    const auto original = std::make_pair(x.pack(), w.pack());
    auto x_inv = x.inverse();
    auto w_inv = w.inverse();
    const auto inverse = std::make_pair(x_inv.pack(), w_inv.pack());
    if (std::make_pair(original.second, original.first) >
        std::make_pair(inverse.second, inverse.first)) {
        x = x_inv;
        w = w_inv;
        return inverse;
    }
    return original;
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

    const auto key = canonicalize_pair(x, w);
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

    const auto key = canonicalize_pair(x, w);
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

    long long result = klc::detail::checked_add(p_x_ws.get_coeff((len_diff - 3) / 2), mu_xs_ws);
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
                current_mu = klc::detail::checked_subtract(current_mu, mu_x_z);
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
    current_mu = klc::detail::checked_subtract(current_mu,
        klc::detail::checked_multiply(mu_z_ws, mu_x_z));
}
