#include "Utils.hpp"
#include <algorithm>

/**
 * @brief Retrieves the pattern formed by a given vector of integers.
 * 
 * Normalizes the elements of the vector into a standard permutation pattern
 * using the numbers 0 to N-1 while preserving their relative order.
 * 
 * @param vec A vector of integers.
 * @return A Permutation object representing the relative order pattern.
 */
Permutation get_pattern(std::vector<int> vec) {
    if (vec.empty()) return Permutation();
    std::vector<int> sorted = vec;
    std::sort(sorted.begin(), sorted.end());
    std::vector<int> res;
    for (int x : vec) {
        // Find the relative rank of each element
        res.push_back(static_cast<int>(std::lower_bound(sorted.begin(), sorted.end(), x) - sorted.begin()));
    }
    return Permutation(res);
}
