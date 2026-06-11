#ifndef UTILS_HPP
#define UTILS_HPP

#include <vector>
#include <string>
#include <cstddef>
#include "Permutation.hpp"

/**
 * @brief Prints a simple progress bar to the standard output.
 * 
 * @param current The current progress value.
 * @param total The total target value for completion.
 */
void print_progress(long long current, long long total);

/**
 * @brief Calculates the factorial of a given integer.
 * 
 * @param n The integer to calculate the factorial for.
 * @return The factorial of n.
 */
long long factorial(int n);

/**
 * @brief A custom hasher structure for the Permutation class.
 * 
 * This enables Permutation objects to be used as keys in standard unordered
 * containers like std::unordered_map or std::unordered_set.
 */
struct PermutationHasher {
    /**
     * @brief Computes the hash value for a given permutation.
     * 
     * @param p The permutation to hash.
     * @return The resulting hash value.
     */
    std::size_t operator()(const Permutation& p) const;
};

/**
 * @brief Retrieves the pattern formed by a given vector of integers.
 * 
 * This normalizes the elements of the vector into a standard permutation pattern
 * using the numbers 0 to N-1 while preserving their relative order.
 * 
 * @param vec A vector of integers.
 * @return A Permutation object representing the relative order pattern.
 */
Permutation get_pattern(std::vector<int> vec);

#endif // UTILS_HPP
