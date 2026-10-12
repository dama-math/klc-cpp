#include "Cli.hpp"
#include "Commands.hpp"
#include <charconv>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
void print_help(const char* program) {
        std::cout << "Usage: " << program << " [options] \n"
                  << "Options:\n"
                  << "  --speedtest                       Speed test KL polynomial calculations\n"
                  << "  --speedtest2 [n] [num] [seed]     Speed test with random permutations (default: n=12, num=1000, seed=42)\n"
                  << "  --random n [num] [seed]           Perform random KL calculations for S_n\n"
                  << "  --verify filename                 Verify KL polynomial calculations against a file\n"
                  << "  --cyclic n                        Run a fast cyclic permutation test for S_n\n"
                  << "  --file filename                   Process permutations from a file\n"
                  << "  --simple w_permutation            Calculate KL polynomial P_e_w(1) for given permutation w\n"
                  << "  x_permutation w_permutation       Calculate KL polynomial for given permutations\n";
}

int parse_integer(std::string_view text) {
    // Preserve acceptance of a leading '+', but require the entire argument.
    if (!text.empty() && text.front() == '+') text.remove_prefix(1);
    int value = 0;
    auto [end, error] = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || error != std::errc{} || end != text.data() + text.size())
        throw std::invalid_argument("Invalid integer argument: " + std::string(text));
    return value;
}

void require_arguments(int argc, int minimum, int maximum, const char* program,
                       std::string_view usage) {
    if (argc < minimum || argc > maximum)
        throw std::invalid_argument("Usage: " + std::string(program) + " " + std::string(usage));
}

int dispatch(int argc, char* argv[]) {
    const std::string_view command = argc > 1 ? argv[1] : "";
    if (command == "--help") {
        require_arguments(argc, 2, 2, argv[0], "--help");
        print_help(argv[0]);
    } else if (command == "--speedtest") {
        require_arguments(argc, 2, 2, argv[0], "--speedtest");
        speedtest();
    } else if (command == "--speedtest2" || command == "--random") {
        const bool benchmark = command == "--speedtest2";
        require_arguments(argc, benchmark ? 2 : 3, 5, argv[0],
            benchmark ? "--speedtest2 [n] [num_calculations] [seed]" :
                        "--random n [num_calculations] [seed]");
        const int n = argc > 2 ? parse_integer(argv[2]) : 12;
        const int count = argc > 3 ? parse_integer(argv[3]) : (benchmark ? 1000 : 10);
        const int seed = argc > 4 ? parse_integer(argv[4]) :
            (benchmark ? 42 : static_cast<int>(std::random_device{}()));
        if (benchmark) speedtest2(n, count, seed);
        else random_kl_calculations(n, count, seed);
    } else if (command == "--verify" || command == "--file") {
        require_arguments(argc, 3, 3, argv[0],
            command == "--verify" ? "--verify filename" : "--file filename");
        KLPolynomialCalculator calculator;
        const bool success = command == "--verify" ? verify_from_file(argv[2], calculator) :
                                                    process_file(argv[2], calculator);
        return success ? 0 : 1;
    } else if (command == "--cyclic") {
        require_arguments(argc, 3, 3, argv[0], "--cyclic n");
        cyclic_test(parse_integer(argv[2]));
    } else if (command == "--simple") {
        require_arguments(argc, 3, 3, argv[0], "--simple w_permutation");
        const Permutation w(argv[2]);
        KLPolynomialCalculator calculator;
        const auto p = calculator.calculate_kl(Permutation::identity(w.size()), w);
        std::cout << p.substitute(1) << std::endl;
    } else {
        require_arguments(argc, 3, 3, argv[0],
                          "x_permutation w_permutation (or --help for options)");
        if (command.starts_with("--"))
            throw std::invalid_argument("Unknown option: " + std::string(command));
        const Permutation x(argv[1]), w(argv[2]);
        KLPolynomialCalculator calculator;
        const auto p = calculator.calculate_kl(x, w);
        std::cout << x.to_string() << " " << w.to_string() << " "
                  << p.to_string_simple() << std::endl;
    }
    return 0;
}
} // namespace

int run_command_line(int argc, char* argv[]) {
    try {
        return dispatch(argc, argv);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
}
