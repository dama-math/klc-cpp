#ifndef COMMANDS_HPP
#define COMMANDS_HPP

#include <string>
#include "KLCalculator.hpp"

/**
 * @brief Process permutations from a file and calculate KL polynomials.
 * @param filename Name of the file containing x and w permutations.
 * @param calculator KLPolynomialCalculator instance.
 */
void process_file(const std::string& filename, KLPolynomialCalculator& calculator);

/**
 * @brief Run a speed test for KL polynomial calculations.
 */
void speedtest();

/**
 * @brief Perform random KL calculations in S_n.
 * @param n Size of the symmetric group.
 * @param num_calculations Number of random calculations to perform.
 * @param seed Random seed.
 */
void random_kl_calculations(int n, int num_calculations, int seed);

/**
 * @brief Verify KL polynomial calculations against a file.
 * @param filename Name of the file containing x, w, and expected polynomials.
 * @param calculator KLPolynomialCalculator instance.
 */
void verify_from_file(const std::string& filename, KLPolynomialCalculator& calculator);

/**
 * @brief Run a fast cyclic permutation test for a specific n.
 * @param n Size of the symmetric group.
 */
void cyclic_test(int n);

#endif // COMMANDS_HPP
