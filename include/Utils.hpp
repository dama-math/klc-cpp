#pragma once

#include <vector>
#include "Permutation.hpp"



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

