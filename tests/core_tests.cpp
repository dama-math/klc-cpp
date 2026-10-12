#include "Bruhat.hpp"
#include "KLCalculator.hpp"
#include "Lascoux.hpp"
#include "detail/FlatMap.hpp"
#include "detail/IntervalSearch.hpp"
#include <algorithm>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
void require(bool condition, const char* message) {
    if (!condition) throw std::runtime_error(message);
}

template<class Exception, class Function>
void throws(Function function) {
    try { function(); }
    catch (const Exception&) { return; }
    throw std::runtime_error("Expected exception was not thrown");
}

std::vector<Permutation> permutations(int n) {
    std::vector<int> data(n);
    std::iota(data.begin(), data.end(), 0);
    std::vector<Permutation> result;
    do { result.emplace_back(data); }
    while (std::next_permutation(data.begin(), data.end()));
    return result;
}

// Independent definition: count entries <= threshold in every prefix.
bool reference_bruhat(const Permutation& x, const Permutation& w) {
    if (x.size() != w.size()) return false;
    for (int prefix = 0; prefix < x.size(); ++prefix)
        for (int threshold = 0; threshold < x.size(); ++threshold) {
            int a = 0, b = 0;
            for (int i = 0; i <= prefix; ++i) {
                a += x[i] <= threshold;
                b += w[i] <= threshold;
            }
            if (a < b) return false;
        }
    return true;
}

// Simple recurrence oracle, without pair reduction, inverse canonicalization or
// extremal pruning. Exhaustively enumerate the correction terms in S_5.
class ReferenceKL {
    std::vector<Permutation> all = permutations(5);
    std::map<std::pair<std::string, std::string>, Polynomial> cache;
public:
    Polynomial calculate(const Permutation& x, const Permutation& w) {
        if (!reference_bruhat(x, w)) return {};
        if (x == w) return Polynomial({1});
        const auto key = std::make_pair(x.to_string(), w.to_string());
        if (auto entry = cache.find(key); entry != cache.end()) return entry->second;
        int s = 0;
        while (w[s] < w[s + 1]) ++s;
        const auto ws = w.multiply_simple_transposition(s);
        const auto xs = x.multiply_simple_transposition(s);
        const int c = xs.length() < x.length() ? 1 : 0;
        Polynomial result;
        result.add_shifted(calculate(xs, ws), 1 - c);
        result.add_shifted(calculate(x, ws), c);
        for (const auto& z : all) {
            const int difference = ws.length() - z.length();
            if (difference <= 0 || difference % 2 == 0 || z[s] < z[s + 1] ||
                !reference_bruhat(x, z) || !reference_bruhat(z, ws)) continue;
            const auto mu = calculate(z, ws).get_coeff((difference - 1) / 2);
            if (mu) result.subtract_shifted(calculate(x, z), (difference + 1) / 2, mu);
        }
        return cache.emplace(key, result).first->second;
    }
};

void test_values() {
    throws<std::invalid_argument>([] { Permutation p(std::vector<int>{0, 0}); });
    throws<std::invalid_argument>([] { Permutation p(std::vector<int>{0, 257}); });
    throws<std::invalid_argument>([] { Permutation p(std::vector<int>{-1, 0}); });
    throws<std::invalid_argument>([] { Permutation p("1234"); });
    throws<std::invalid_argument>([] { Permutation p(nullptr, 2); });
    throws<std::length_error>([] { Permutation p(nullptr, -1); });
    throws<std::length_error>([] { Permutation::identity(Permutation::MAX_SIZE + 1); });
    throws<std::out_of_range>([] { Permutation("01").multiply_simple_transposition(1); });
    throws<std::out_of_range>([] { Permutation("01").at(-1); });
    throws<std::out_of_range>([] { Permutation("01").at(2); });
    require(Permutation(nullptr, 0).size() == 0, "Empty permutation");
    auto p = Polynomial({1, 2});
    p.add_shifted(p, 1);
    require(p == Polynomial({1, 3, 2}), "Aliased polynomial addition");
    p.subtract_shifted(p, 0);
    require(p.is_zero(), "Aliased polynomial subtraction");
    throws<std::invalid_argument>([] { Polynomial({1}).add_shifted(Polynomial({1}), -1); });
    throws<std::length_error>([] { Polynomial({1}).add_shifted(Polynomial({1}), Polynomial::MAX_CAPACITY); });
    constexpr auto maximum = std::numeric_limits<long long>::max();
    constexpr auto minimum = std::numeric_limits<long long>::min();
    throws<std::overflow_error>([=] { Polynomial({maximum}).add_shifted(Polynomial({1}), 0); });
    throws<std::overflow_error>([=] { Polynomial({minimum}).subtract_shifted(Polynomial({1}), 0); });
    throws<std::overflow_error>([=] { Polynomial p({maximum}); p *= Polynomial({2}); });
    throws<std::overflow_error>([=] { Polynomial({maximum, 1}).substitute(1); });
    require(Polynomial({1}).substitute(maximum) == 1, "Constant evaluation must not overflow");
    require(Polynomial({1, 2, 3}).substitute(2) == 17, "Polynomial evaluation");
    require(klc::detail::checked_multiply(minimum, 1) == minimum, "Signed multiplication boundary");
    throws<std::overflow_error>([=] { klc::detail::checked_multiply(minimum, -1); });
    KLPolynomialCalculator calculator;
    throws<std::invalid_argument>([&] { calculator.calculate_kl(Permutation("0"), Permutation("01")); });
    throws<std::invalid_argument>([] { auto x = Permutation("0"), w = Permutation("01"); Bruhat::reduce(x, w); });
    throws<std::invalid_argument>([] { Lascoux::calculate_kl(Permutation("0"), Permutation("01")); });
    throws<std::invalid_argument>([] { Lascoux::calculate_kl(Permutation("0123"), Permutation("2301")); });
}

void test_cache() {
    struct CollisionHasher { size_t operator()(int) const { return 0; } };
    klc::detail::FlatMap<int, int, CollisionHasher> map(1);
    for (int i = 0; i < 200; ++i) map.insert_or_assign(i, i * 3);
    for (int i = 0; i < 200; ++i) require(*map.find(i) == i * 3, "Cache growth/collision");
    map.insert_or_assign(10, 42);
    require(*map.find(10) == 42 && map.size() == 200, "Cache update");
    require(map.find(201) == nullptr, "Missing cache key");
    map.clear();
    require(map.find(10) == nullptr && map.size() == 0, "Cache clear");
    klc::detail::FlatMap<int, Polynomial, CollisionHasher> aliased(1);
    aliased.insert_or_assign(0, Polynomial({1, 2, 3}));
    aliased.insert_or_assign(1, *aliased.find(0));
    require(*aliased.find(1) == Polynomial({1, 2, 3}), "Aliased cache insertion during growth");
}

void test_algorithms() {
    const auto all = permutations(5);
    ReferenceKL reference;
    KLPolynomialCalculator calculator;
    for (const auto& x : all)
        for (const auto& w : all) {
            const bool expected_order = reference_bruhat(x, w);
            require(Bruhat::is_bruhat_less_than(x, w) == expected_order, "Bruhat oracle mismatch");
            require(klc::detail::scalar_bruhat_compare(x, w) == expected_order, "Scalar Bruhat mismatch");
            const auto expected = reference.calculate(x, w);
            const auto result = calculator.calculate_kl(x, w, false);
            require(result == expected, "KL recurrence oracle mismatch");
            require(calculator.calculate_kl(x.inverse(), w.inverse(), false) == result, "Inverse symmetry");
            if (expected_order) {
                auto reduced_x = x, reduced_w = w;
                Bruhat::reduce(reduced_x, reduced_w);
                require(calculator.calculate_kl(reduced_x, reduced_w, false) == result, "Reduction invariance");
            }
        }
    calculator.clear_cache();
    require(calculator.calculate_kl(Permutation("0123"), Permutation("2301")) == Polynomial({1, 1}),
            "Known KL value after cache clear");

    std::mt19937 rng(42);
    for (int n : {0, 1, 8, 16, Permutation::MAX_SIZE}) {
        auto x_data = Permutation::identity(n).to_vector(), w_data = x_data;
        for (int trial = 0; trial < 100; ++trial) {
            std::shuffle(x_data.begin(), x_data.end(), rng);
            std::shuffle(w_data.begin(), w_data.end(), rng);
            Permutation x(x_data), w(w_data);
            require(Bruhat::is_bruhat_less_than(x, w) == reference_bruhat(x, w), "Boundary Bruhat comparison");
            require(x.inverse().inverse() == x, "Permutation inverse");
        }
    }
    const int n = Permutation::MAX_SIZE;
    auto x_data = Permutation::identity(n).to_vector();
    auto w_data = x_data;
    w_data[0] = 2; w_data[1] = 3; w_data[2] = 0; w_data[3] = 1;
    require(calculator.calculate_kl(Permutation(x_data), Permutation(w_data), false) == Polynomial({1, 1}),
            "Maximum-size embedded KL calculation");
    for (int trial = 0; trial < 300; ++trial) {
        auto a = Permutation::identity(8).to_vector(), b = a;
        std::shuffle(a.begin(), a.end(), rng);
        std::shuffle(b.begin(), b.end(), rng);
        const Permutation x(a), w(b);
        if (!Lascoux::has_3412_pattern(w))
            require(calculator.calculate_kl(x, w) == calculator.calculate_kl(x, w, false),
                    "Lascoux/recurrence agreement");
    }
}

void test_search() {
    auto all = permutations(5);
    for (size_t trial = 0; trial < 120; ++trial) {
        const auto& x = all[trial];
        const auto& ws = all[(trial * 37) % all.size()];
        const int s = static_cast<int>(trial % 4);
        std::set<std::string> actual, expected;
        klc::detail::search_extremal_pairs(x, ws, s, [&](const Permutation& z) {
            actual.insert(z.to_string()); return true;
        });
        for (const auto& z : all)
            if (reference_bruhat(x, z) && reference_bruhat(z, ws) &&
                z[s] > z[s + 1] && Bruhat::is_extremal_pair(z, ws)) expected.insert(z.to_string());
        require(actual == expected, "Interval enumeration oracle mismatch");
        std::set<std::string> scalar;
        klc::detail::enumerate_extremal_pairs<klc::detail::ScalarRank>(x, ws, s,
            [&](const Permutation& z) { scalar.insert(z.to_string()); return true; });
        require(scalar == expected, "Scalar interval enumeration mismatch");
        int calls = 0;
        klc::detail::search_extremal_pairs(x, ws, s, [&](const Permutation&) { ++calls; return false; });
        require(calls == (expected.empty() ? 0 : 1), "Search early termination");
    }
    int calls = 0;
    auto callback = [&](const Permutation&) { ++calls; return true; };
    klc::detail::search_extremal_pairs(Permutation(), Permutation(), -1, callback);
    klc::detail::search_extremal_pairs(Permutation("0"), Permutation("0"), -1, callback);
    require(calls == 1, "Small-size search");
    // Exercise the >16 backend without enumerating a huge interval.
    auto identity = Permutation::identity(Permutation::MAX_SIZE);
    calls = 0;
    klc::detail::search_extremal_pairs(identity, identity, 0, callback);
    require(calls == 0, "Maximum-size search");
    auto adjacent = identity.multiply_simple_transposition(identity.size() - 2);
    calls = 0;
    klc::detail::search_extremal_pairs(identity, adjacent, identity.size() - 2,
        [&](const Permutation& z) {
            require(z == adjacent && z.length() == 1, "Maximum-size candidate");
            ++calls; return true;
        });
    require(calls == 1, "Maximum-size nonempty search");
}
} // namespace

int main() {
    try {
        test_values(); test_cache(); test_algorithms(); test_search();
        std::cout << "Core tests passed (S_5: 14,400 pairs, independent recurrence and rank oracles).\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Test failed: " << e.what() << '\n';
        return 1;
    }
}
