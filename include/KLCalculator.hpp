#pragma once

#include "Config.hpp"
#include "Permutation.hpp"
#include "Polynomial.hpp"
#include "detail/FlatMap.hpp"
#include <utility>


/**
 * @brief High-performance Kazhdan-Lusztig polynomial calculator.
 * 
 * Computes KL polynomials P_{x, w}(q) using Warrington's recursive algorithm
 * with pair reductions, SIMD-accelerated extremal pair searches, covering relation
 * filtering, and memoization. Also supports Lascoux's fast algorithm for vexillary
 * (3412-avoiding) permutations.
 */
class KLPolynomialCalculator {
public:
    KLPolynomialCalculator();

    /**
     * @brief Computes the Kazhdan-Lusztig polynomial P_{x, w}(q).
     * 
     * @param x Inner permutation.
     * @param w Outer permutation.
     * @param use_lascoux If true, uses Lascoux's method when w is 3412-avoiding.
     * @return The computed polynomial (0 if x is not <= w in Bruhat order).
     */
    Polynomial calculate_kl(const Permutation& x, const Permutation& w, bool use_lascoux = true);

    /** @brief Clears memoization tables for polynomials and mu coefficients. */
    void clear_cache();

private:
    struct FastPermPairHasher {
        static inline uint64_t splitmix64(uint64_t x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        std::size_t operator()(const std::pair<PackedInt, PackedInt>& p) const {
#ifdef USE_128BIT_PACKING
            uint64_t h1 = static_cast<uint64_t>(p.first ^ (p.first >> 64));
            uint64_t h2 = static_cast<uint64_t>(p.second ^ (p.second >> 64));
#else
            uint64_t h1 = p.first;
            uint64_t h2 = p.second;
#endif
            return splitmix64(h1 ^ splitmix64(h2));
        }
    };

    using PermPair = std::pair<PackedInt, PackedInt>;

    // Memoization table for full polynomials P_{x, w}(q)
    klc::detail::FlatMap<PermPair, Polynomial, FastPermPairHasher> cache;
    
    // Memoization table for leading mu coefficients mu(x, w)
    klc::detail::FlatMap<PermPair, long long, FastPermPairHasher> mu_cache;

    // --- Core Recurrences ---
    Polynomial get_kl_recursive(Permutation x, Permutation w, bool already_bruhat_checked);
    long long get_mu_recursive(Permutation x, Permutation w, bool already_bruhat_checked);

    // --- Correction Terms for P_{x, w}(q) ---
    void apply_correction_term(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void find_extremal_pairs_and_correct(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void find_covering_relations_and_correct(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void perform_correction(Polynomial& poly, const Permutation& x, const Permutation& ws, const Permutation& z);

    // --- Correction Terms for mu(x, w) ---
    void mu_apply_correction_term(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_find_extremal_pairs_and_correct(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_find_covering_relations_and_correct(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_perform_correction(long long& current_mu, const Permutation& x, const Permutation& ws, const Permutation& z);
};
