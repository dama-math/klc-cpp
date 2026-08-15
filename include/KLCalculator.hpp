#pragma once

#include "Config.hpp"
#include "Permutation.hpp"
#include "Polynomial.hpp"
#include <unordered_map>
#include <string>
#include <vector>
#include <utility> // For std::pair

/**
 * @brief High-performance Kazhdan-Lusztig polynomial calculator.
 * 
 * Implements various algorithms including the standard recursive formula and
 * Lascoux's optimization for pattern-free permutations. Employs advanced
 * caching (interning, fast maps) to eliminate redundant calculations.
 */
class KLPolynomialCalculator {
public:
    /** @brief Default constructor. */
    KLPolynomialCalculator();

    /**
     * @brief Calculates the Kazhdan-Lusztig polynomial P_{x, w}(q).
     * 
     * @param x Inner permutation.
     * @param w Outer permutation.
     * @param use_lascoux If true, attempt to use Lascoux's algorithm for pattern-free permutations.
     * @return The computed Kazhdan-Lusztig polynomial.
     */
    Polynomial calculate_kl(const Permutation& x, const Permutation& w, bool use_lascoux = true);

    /** 
     * @brief Clears the internal computation caches.
     * 
     * Frees all stored polynomials and partially computed values to reclaim memory.
     */
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

    // Fast Cache mapping packed permutation pairs to polynomials
    std::unordered_map<PermPair, Polynomial, FastPermPairHasher> cache; 
    
    // Fast Cache for tracking intermediate mu coefficients
    std::unordered_map<PermPair, long long, FastPermPairHasher> mu_cache;

    /** @brief Core recursive function for standard KL calculation. */
    Polynomial get_kl_recursive(Permutation x, Permutation w, bool already_bruhat_checked);

    /** @brief Coordinates the application of all correction terms for a recursive step. */
    void apply_correction_term(
        Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Finds and applies correction polynomials derived from extremal pairs. */
    void find_extremal_pairs_and_correct(
        Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Finds and applies correction polynomials derived from covering relations. */
    void find_covering_relations_and_correct(
        Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Performs polynomial subtraction for a specific valid element z. */
    void perform_correction(
        Polynomial& poly, const Permutation& x, const Permutation& ws, const Permutation& z);

    /** @brief Retrieves the raw mu coefficient for a specific pair. */
    long long get_mu_coefficient(const Permutation& z, const Permutation& ws);

    /** @brief Core recursive function dedicated to computing mu coefficients rapidly. */
    long long get_mu_recursive(Permutation x, Permutation w, bool already_bruhat_checked);

    /** @brief Orchestrates corrections during mu coefficient derivation. */
    void mu_apply_correction_term(
        long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Corrects mu based on relevant extremal pairs. */
    void mu_find_extremal_pairs_and_correct(
        long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Corrects mu based on relevant covering relations. */
    void mu_find_covering_relations_and_correct(
        long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);

    /** @brief Subtracts the mu product correction for a specific valid element z. */
    void mu_perform_correction(
        long long& current_mu, const Permutation& x, const Permutation& ws, const Permutation& z);

};
