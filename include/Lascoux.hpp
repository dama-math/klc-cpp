#pragma once

#include "Permutation.hpp"
#include "Polynomial.hpp"

/**
 * @brief Namespace for Lascoux's algorithm and related optimizations.
 */
namespace Lascoux {

/**
 * @brief Checks if a permutation contains a 3412 pattern.
 * 
 * A permutation has a 3412 pattern if there exist indices i < j < k < l
 * such that w[k] < w[l] < w[i] < w[j].
 * 
 * @param p The permutation to check.
 * @return True if the pattern exists, false otherwise.
 */
bool has_3412_pattern(const Permutation& p);

/**
 * @brief Calculates the Kazhdan-Lusztig polynomial using Lascoux's algorithm.
 * 
 * This algorithm is specifically optimized for permutations that do not contain
 * a 3412 pattern. It significantly speeds up computations for such permutations.
 * 
 * @param x The inner permutation.
 * @param w The outer permutation.
 * @return The calculated Kazhdan-Lusztig polynomial.
 */
Polynomial calculate_kl(const Permutation& x, const Permutation& w);

} // namespace Lascoux
