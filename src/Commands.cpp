#include "Commands.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <algorithm>
#include <chrono>
#include <random>
#include <iomanip>

void process_file(const std::string& filename, KLPolynomialCalculator& calculator) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }

    std::string line;
    while (std::getline(file, line)) {
        size_t first_non_space = line.find_first_not_of(" \t\r\n");
        if (first_non_space == std::string::npos || line[first_non_space] == '#') continue;

        std::stringstream ss(line);
        std::string x_str, w_str;
        if (!(ss >> x_str >> w_str)) continue;

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

void speedtest() {
    constexpr int MAX_N_SPEED = 9;
    std::cout << "Calculating KLs for n <= " << MAX_N_SPEED << "..." << std::endl;
    KLPolynomialCalculator calculator;

    for (int n = 7; n <= MAX_N_SPEED; ++n) {
        std::vector<int> p_data(n);
        std::iota(p_data.begin(), p_data.end(), 0);
        Permutation e(p_data);

        auto start = std::chrono::high_resolution_clock::now();
        
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

void speedtest2(int n, int num_calculations, int seed) {
    std::cout << "Speed test 2: Calculating P_{e,w}(q) for " << num_calculations
              << " random permutations in S_" << n << " (seed: " << seed << ")..." << std::endl;

    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    Permutation e(p);

    std::mt19937 rng(seed);
    std::vector<Permutation> permutations;
    permutations.reserve(num_calculations);
    for (int i = 0; i < num_calculations; ++i) {
        std::shuffle(p.begin(), p.end(), rng);
        permutations.emplace_back(p);
    }

    KLPolynomialCalculator calculator;
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < num_calculations; ++i) {
        calculator.calculate_kl(e, permutations[i]);
    }
    auto end = std::chrono::high_resolution_clock::now();

    std::chrono::duration<double> diff = end - start;
    std::cout << "Calculation finished." << std::endl;
    std::cout << "Total time: " << std::fixed << std::setprecision(4) << diff.count() << " s" << std::endl;
    if (num_calculations > 0) {
        std::cout << "Average time: " << std::fixed << std::setprecision(4)
                  << (diff.count() / num_calculations) * 1000.0 << " ms/calc" << std::endl;
    }
}

void random_kl_calculations(int n, int num_calculations, int seed) {
    KLPolynomialCalculator calculator;
    
    std::vector<int> p(n);
    std::iota(p.begin(), p.end(), 0);
    Permutation e(p);

    std::mt19937 rng(seed);
    std::cout << "Performing " << num_calculations << " random KL calculations in S_" << n 
              << " with seed " << seed << "..." << std::endl;

    for (int i = 0; i < num_calculations; ++i) {
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
        size_t first_non_space = line.find_first_not_of(" \t\r\n");
        if (first_non_space == std::string::npos || line[first_non_space] == '#') continue;

        std::stringstream ss(line);
        std::string x_str, w_str, expected_poly;
        if (!(ss >> x_str >> w_str >> expected_poly)) continue;

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
