#pragma once

#include <string>
#include "KLCalculator.hpp"

/**
 * @brief Process permutations from a file and calculate KL polynomials.
 * @param filename Name of the file containing x and w permutations.
 * @param calculator KLPolynomialCalculator instance.
 * @return False if a line could not be processed or the file could not be opened.
 */
bool process_file(const std::string& filename, KLPolynomialCalculator& calculator);

/**
 * @brief Run a speed test for KL polynomial calculations.
 */
void speedtest();

/**
 * @brief Run a speed test calculating P_{e, w}(q) for random permutations.
 * @param n Size of the symmetric group (default: 12).
 * @param num_calculations Number of random calculations to perform (default: 1000).
 * @param seed Random seed (default: 42).
 */
void speedtest2(int n = 12, int num_calculations = 1000, int seed = 42);

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
 * @return False on a mismatch, malformed input or file-open failure.
 */
bool verify_from_file(const std::string& filename, KLPolynomialCalculator& calculator);

/**
 * @brief Run a fast cyclic permutation test for a specific n.
 * @param n Size of the symmetric group.
 */
void cyclic_test(int n);

