#include "Commands.hpp"
#include "Utils.hpp"
#include "Lascoux.hpp"
#include <iostream>
#include <fstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include <iomanip>

/**
 * @brief Processes a file containing pairs of permutations and calculates their KL polynomials.
 * 
 * Each line in the file should contain two space-separated permutations, x and w.
 * The function reads the file line by line, computes the KL polynomial P_{x,w}, 
 * and outputs the result to the standard output.
 * 
 * @param filename The path to the file to be processed.
 * @param calculator The KL polynomial calculator instance to use.
 */
void process_file(const std::string& filename, KLPolynomialCalculator& calculator) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        // Find the space separating the two permutations
        size_t space_pos = line.find(' ');
        if (space_pos == std::string::npos) continue;

        std::string x_str = line.substr(0, space_pos);
        std::string w_str = line.substr(space_pos + 1);

        try {
            Permutation x(x_str);
            Permutation w(w_str);

            std::cout << x.to_string() << " " << w.to_string() << " ";
            Polynomial p = calculator.calculate_kl(x, w);
            std::cout << p.to_string_simple() << std::endl;

        } catch (const std::exception& e) {
            std::cerr << "Error processing line: " << line << " -> " << e.what() << std::endl;
        }
    }
}

/**
 * @brief Runs a speed test for KL polynomial calculations.
 * 
 * This function calculates the KL polynomial for all pairs (e, w) in the 
 * symmetric groups S_7 up to S_9, where e is the identity permutation.
 * It outputs the time taken for each group.
 */
void speedtest() {
    const int MAX_N_SPEED = 9;
    std::cout << "Calculating KLs for n <= " << MAX_N_SPEED << "..." << std::endl;
    KLPolynomialCalculator calculator;

    for (int n = 7; n <= MAX_N_SPEED; ++n) {
        // Create the identity permutation
        std::vector<int> p_data(n);
        std::iota(p_data.begin(), p_data.end(), 0);
        Permutation e(p_data);

        auto start = std::chrono::high_resolution_clock::now();
        
        // Iterate through all permutations in S_n
        do {
            Permutation w(p_data);
            calculator.calculate_kl(e, w);
        } while (std::next_permutation(p_data.begin(), p_data.end()));
        
        auto end = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = end - start;
        std::cout << "Calculation for n=" << n << " finished. Time: " << diff.count() << " s" << std::endl;
    }
    std::cout << "Calculation finished." << std::endl;
}

/**
 * @brief Performs a specified number of random KL polynomial calculations.
 * 
 * This function generates random permutations in S_n and computes 
 * their KL polynomials relative to the identity permutation.
 * 
 * @param n The size of the symmetric group.
 * @param num_calculations The number of random calculations to perform.
 * @param seed The random seed for reproducibility.
 */
void random_kl_calculations(int n, int num_calculations, int seed) {
    KLPolynomialCalculator calculator;
    
    // Create the identity permutation
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    Permutation e(p);

    std::mt19937 rng(seed);
    std::cout << "Performing " << num_calculations << " random KL calculations in S_" << n 
              << " with seed " << seed << "..." << std::endl;

    for (int i = 0; i < num_calculations; ++i) {
        // Generate a random permutation
        std::shuffle(p.begin(), p.end(), rng);
        Permutation w(p);
        
        int length = w.length();
        Polynomial P_e_w = calculator.calculate_kl(e, w);
        
        std::cout << "Calc " << std::setw(3) << (i + 1) 
              << ": w = " << std::setw(n) << std::left << w.to_string()
              << ", l(w) = " << std::setw(3) << std::right << length
              << ", P_e_w(q) = " << P_e_w.to_string_simple() << std::endl;
    }
}

void verify_from_file(const std::string& filename, KLPolynomialCalculator& calculator) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    int passed = 0, total = 0;
    while (std::getline(file, line)) {
        if (line.empty()) continue;
        size_t sp1 = line.find(' ');
        size_t sp2 = line.find(' ', sp1 + 1);
        if (sp1 == std::string::npos || sp2 == std::string::npos) continue;

        std::string x_str = line.substr(0, sp1);
        std::string w_str = line.substr(sp1 + 1, sp2 - sp1 - 1);
        std::string expected_poly = line.substr(sp2 + 1);

        try {
            Permutation x(x_str);
            Permutation w(w_str);
            Polynomial p = calculator.calculate_kl(x, w);
            std::string result_poly = p.to_string_simple();
            
            total++;
            if (result_poly == expected_poly) {
                passed++;
                std::cout << "[PASS] " << x_str << " " << w_str << std::endl;
            } else {
                std::cout << "[FAIL] " << x_str << " " << w_str 
                          << "\n  Expected: " << expected_poly 
                          << "\n  Got:      " << result_poly << std::endl;
            }
        } catch (const std::exception& e) {
            std::cerr << "Error processing line: " << line << " -> " << e.what() << std::endl;
        }
    }
    std::cout << "Verification complete: " << passed << "/" << total << " passed." << std::endl;
}

void cyclic_test(int n) {
    if (n < 2) return;
    KLPolynomialCalculator calculator;
    
    std::vector<int> e_data(n);
    std::iota(e_data.begin(), e_data.end(), 0);
    Permutation e(e_data);

    std::vector<int> w_data(n);
    for (int i = 0; i < n - 1; ++i) {
        w_data[i] = i + 1;
    }
    w_data[n - 1] = 0;
    Permutation w(w_data);

    std::cout << "Performing cyclic test for n=" << n << std::endl;
    std::cout << "w = " << w.to_string() << std::endl;

    auto start = std::chrono::high_resolution_clock::now();
    Polynomial p = calculator.calculate_kl(e, w);
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Result: " << p.to_string_simple() << std::endl;
    std::cout << "Time: " << diff.count() << " s" << std::endl;
}
