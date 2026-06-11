#include "Utils.hpp"
#include <iostream>
#include <iomanip>
#include <numeric>
#include <algorithm>
#include <functional>

/**
 * @brief Prints a simple progress bar to the standard output.
 * 
 * @param current The current progress value.
 * @param total The total target value for completion.
 */
void print_progress(long long current, long long total) {
    int width = 50;
    double progress = static_cast<double>(current) / total;
    int pos = static_cast<int>(width * progress);

    std::cout << "\r[";
    for (int i = 0; i < width; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << std::fixed << std::setprecision(1) << progress * 100.0 << "% ("
              << current << "/" << total << ")" << std::flush;
}

/**
 * @brief Calculates the factorial of a given integer.
 * 
 * @param n The integer to calculate the factorial for.
 * @return The factorial of n.
 */
long long factorial(int n) {
    long long res = 1;
    for (int i = 2; i <= n; ++i) res *= i;
    return res;
}

/**
 * @brief Computes the hash value for a given permutation.
 * 
 * @param p The permutation to hash.
 * @return The resulting hash value.
 */
std::size_t PermutationHasher::operator()(const Permutation& p) const {
    std::size_t seed = p.size();
    for (auto val : p) {
        // val is int8_t, cast to int for stable hashing
        seed ^= std::hash<int>()((int)val) + 0x9e3779b9 + (seed << 6) + (seed >> 8);
    }
    return seed;
}

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
