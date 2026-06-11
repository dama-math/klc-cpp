#pragma once

#include "Permutation.hpp"

/**
 * @brief Namespace for operations related to the Bruhat order on permutations.
 */
namespace Bruhat {

/**
 * @brief Reduces the permutation pair (x, w) in place.
 * 
 * This reduction applies standard transformations to simplify the pair 
 * without changing their Bruhat order relationship or associated KL polynomials.
 * Corresponds to olp_reduce_inplace in the original C algorithms.
 * 
 * @param x The first permutation (will be modified in place).
 * @param w The second permutation (will be modified in place).
 */
void reduce(Permutation& x, Permutation& w);

/**
 * @brief Checks if the permutations (x, w) form an extremal pair.
 * 
 * Corresponds to olp_isep in the original C algorithms.
 * 
 * @param x The first permutation.
 * @param w The second permutation.
 * @return True if (x, w) is an extremal pair, false otherwise.
 */
bool is_extremal_pair(const Permutation& x, const Permutation& w);

/**
 * @brief Overload for checking extremal pairs using raw vectors.
 * 
 * @param x_data The first permutation vector.
 * @param w_data The second permutation vector.
 * @return True if (x, w) is an extremal pair, false otherwise.
 */
bool is_extremal_pair(const std::vector<int>& x_data, const std::vector<int>& w_data);

/**
 * @brief Checks if permutation x is less than or equal to permutation w in the Bruhat order.
 * 
 * This is a high-performance implementation utilizing the rank matrix algorithm
 * optimized for smaller permutations.
 * 
 * @param x The first permutation.
 * @param w The second permutation.
 * @return True if x <= w in the Bruhat order, false otherwise.
 */
bool is_bruhat_less_than(const Permutation& x, const Permutation& w);

/**
 * @brief Overload for checking Bruhat order using raw vectors.
 * 
 * @param x_data The first permutation vector.
 * @param w_data The second permutation vector.
 * @return True if x <= w in the Bruhat order, false otherwise.
 */
bool is_bruhat_less_than(const std::vector<int>& x_data, const std::vector<int>& w_data);

} // namespace Bruhat
