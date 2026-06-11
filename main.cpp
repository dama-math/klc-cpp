#include <iostream>
#include <string>
#include <vector>
#include <numeric>
#include <random>
#include "Permutation.hpp"
#include "KLCalculator.hpp"
#include "Commands.hpp"
#include "Utils.hpp"

/**
 * @brief Entry point for the KL polynomial calculator.
 * 
 * Provides various options to compute KL polynomials and perform tests.
 */
int main(int argc, char* argv[]) {
    if (argc > 1 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: " << argv[0] << " [options] \n"
                  << "Options:\n"
                  << "  --speedtest                       Speed test KL polynomial calculations\n"
                  << "  --random n [num] [seed]           Perform random KL calculations for S_n\n"
                  << "  --verify filename                 Verify KL polynomial calculations against a file\n"
                  << "  --cyclic n                        Run a fast cyclic permutation test for S_n\n"
                  << "  --file filename                   Process permutations from a file\n"
                  << "  --simple w_permutation            Calculate KL polynomial P_e_w(1) for given permutation w\n"
                  << "  x_permutation w_permutation       Calculate KL polynomial for given permutations\n";
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--speedtest") {
        speedtest();
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--random") {
        if (argc < 3) {
            std::cerr << "Usage: " << argv[0] << " --random n [num_calculations] [seed]" << std::endl;
            return 1;
        }
        try {
            int n = std::stoi(argv[2]);
            int num_calculations = (argc > 3) ? std::stoi(argv[3]) : 10;
            int seed = (argc > 4) ? std::stoi(argv[4]) : std::random_device{}();
            random_kl_calculations(n, num_calculations, seed);
        } catch (const std::exception& e) {
            std::cerr << "Error: invalid argument: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--verify") {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " --verify filename" << std::endl;
            return 1;
        }
        KLPolynomialCalculator calculator;
        verify_from_file(argv[2], calculator);
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--cyclic") {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " --cyclic n" << std::endl;
            return 1;
        }
        try {
            int n = std::stoi(argv[2]);
            cyclic_test(n);
        } catch (const std::exception& e) {
            std::cerr << "Error: invalid argument: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--file") {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " --file filename" << std::endl;
            return 1;
        }
        KLPolynomialCalculator calculator;
        process_file(argv[2], calculator);
        return 0;
    }

    if (argc > 1 && std::string(argv[1]) == "--simple") {
        if (argc != 3) {
            std::cerr << "Usage: " << argv[0] << " --simple w_permutation" << std::endl;
            return 1;
        }
        KLPolynomialCalculator calculator;
        try {
            Permutation w(argv[2]);
            std::vector<int> identity(w.size());
            std::iota(identity.begin(), identity.end(), 0);
            Permutation x(identity);
            Polynomial P = calculator.calculate_kl(x, w);
            std::cout << P.substitute(1) << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    // Default behavior when two permutations are provided
    if (argc == 3) {
        KLPolynomialCalculator calculator;
        try {
            Permutation x(argv[1]);
            Permutation w(argv[2]);
            std::cout << x.to_string() << " " << w.to_string() << " ";
            Polynomial p = calculator.calculate_kl(x, w);
            std::cout << p.to_string_simple() << std::endl;
        } catch (const std::exception& e) {
            std::cerr << "Error: " << e.what() << std::endl;
            return 1;
        }
        return 0;
    }

    std::cerr << "Usage: " << argv[0] << " x_permutation w_permutation (or --help for options)" << std::endl;
    return 1;
}
