#pragma once

#include "Permutation.hpp"
#include "Polynomial.hpp"

/**
 * @brief Lascoux's algorithm for Kazhdan-Lusztig polynomials of 3412-avoiding (vexillary) permutations.
 */
namespace Lascoux {

/**
 * @brief Checks if permutation p contains a 3412 pattern.
 * 
 * A permutation contains 3412 if there exist indices i < j < k < l such that
 * p[k] < p[l] < p[i] < p[j].
 */
bool has_3412_pattern(const Permutation& p);

/**
 * @brief Computes the KL polynomial P_{x, w}(q) using Lascoux's tree method.
 * 
 * Applicable when w is 3412-avoiding (vexillary).
 */
Polynomial calculate_kl(const Permutation& x, const Permutation& w);

} // namespace Lascoux
