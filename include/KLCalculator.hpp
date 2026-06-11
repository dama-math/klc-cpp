#pragma once

#include "Config.hpp"
#include "Permutation.hpp"
#include "Polynomial.hpp"
#include "PolynomialPool.hpp"
#include "FastMap.hpp"
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

    /** 
     * @brief Precomputes KL polynomials for all permutations up to a given size.
     * 
     * Useful for eagerly filling the cache when batch processing.
     * 
     * @param max_n The maximum permutation size for precomputation.
     */
    void precompute_kl_polynomials(int max_n);

private:
    struct CacheKeyHasher {
        std::size_t operator()(const std::pair<PackedInt, PackedInt>& p) const {
#ifdef USE_128BIT_PACKING
             uint64_t h1 = static_cast<uint64_t>(p.first ^ (p.first >> 64));
             uint64_t h2 = static_cast<uint64_t>(p.second ^ (p.second >> 64));
#else
             uint64_t h1 = p.first;
             uint64_t h2 = p.second;
#endif
             return std::hash<uint64_t>{}(h1) ^ (std::hash<uint64_t>{}(h2) << 1);
        }
    };
    
    using PermPair = std::pair<PackedInt, PackedInt>;
    
    // Interning Pool for polynomials to save memory allocations
    PolynomialPool poly_pool;

    // Fast Cache mapping packed permutation pairs to polynomial IDs
    FastMap<PermPair, uint32_t, CacheKeyHasher> cache; 
    
    // Fast Cache for tracking intermediate mu coefficients
    FastMap<PermPair, long long, CacheKeyHasher> mu_cache;

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

    /** @brief Computes rank matrices for legacy Bruhat order comparison. */
    static std::vector<std::vector<int>> compute_rx_matrix(const std::vector<int>& p_data);
};
