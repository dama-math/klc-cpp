#pragma once

#include "Config.hpp"
#include "Permutation.hpp"
#include "Polynomial.hpp"
#include <unordered_map>
#include <utility>



/**
 * @brief High-performance Kazhdan-Lusztig polynomial calculator.
 * 
 * Computes KL polynomials P_{x, w}(q) using Warrington's recursive algorithm
 * with pair reductions, SIMD-accelerated extremal pair searches, covering relation
 * filtering, and memoization. Also supports Lascoux's fast algorithm for vexillary
 * (3412-avoiding) permutations.
 */
class KLPolynomialCalculator {
public:
    KLPolynomialCalculator();

    /**
     * @brief Computes the Kazhdan-Lusztig polynomial P_{x, w}(q).
     * 
     * @param x Inner permutation.
     * @param w Outer permutation.
     * @param use_lascoux If true, uses Lascoux's method when w is 3412-avoiding.
     * @return The computed polynomial (0 if x is not <= w in Bruhat order).
     */
    Polynomial calculate_kl(const Permutation& x, const Permutation& w, bool use_lascoux = true);

    /** @brief Clears memoization tables for polynomials and mu coefficients. */
    void clear_cache();

private:
    struct FastPermPairHasher {
        static inline uint64_t splitmix64(uint64_t x) {
            x += 0x9e3779b97f4a7c15ULL;
            x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
            x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
            return x ^ (x >> 31);
        }

        std::size_t operator()(const std::pair<PackedInt, PackedInt>& p) const {
#ifdef USE_128BIT_PACKING
            uint64_t h1 = static_cast<uint64_t>(p.first ^ (p.first >> 64));
            uint64_t h2 = static_cast<uint64_t>(p.second ^ (p.second >> 64));
#else
            uint64_t h1 = p.first;
            uint64_t h2 = p.second;
#endif
            return splitmix64(h1 ^ splitmix64(h2));
        }
    };

    using PermPair = std::pair<PackedInt, PackedInt>;

template <typename Key, typename Value, typename Hasher>
class FastFlatMap {
public:
    FastFlatMap(size_t initial_cap = 4096) {
        capacity_ = 1;
        while (capacity_ < initial_cap) capacity_ <<= 1;
        keys_.resize(capacity_);
        values_.resize(capacity_);
        occupied_.assign(capacity_, 0);
    }

    void clear() {
        std::fill(occupied_.begin(), occupied_.end(), 0);
        size_ = 0;
    }

    size_t size() const { return size_; }

    Value* find(const Key& key) {
        if (size_ == 0) return nullptr;
        size_t mask = capacity_ - 1;
        size_t idx = hasher_(key) & mask;
        while (occupied_[idx]) {
            if (keys_[idx] == key) {
                return &values_[idx];
            }
            idx = (idx + 1) & mask;
        }
        return nullptr;
    }

    void insert_or_assign(const Key& key, const Value& val) {
        if (size_ * 10 >= capacity_ * 7) {
            grow();
        }
        size_t mask = capacity_ - 1;
        size_t idx = hasher_(key) & mask;
        while (occupied_[idx]) {
            if (keys_[idx] == key) {
                values_[idx] = val;
                return;
            }
            idx = (idx + 1) & mask;
        }
        occupied_[idx] = 1;
        keys_[idx] = key;
        values_[idx] = val;
        size_++;
    }

private:
    void grow() {
        size_t old_cap = capacity_;
        std::vector<Key> old_keys = std::move(keys_);
        std::vector<Value> old_values = std::move(values_);
        std::vector<uint8_t> old_occupied = std::move(occupied_);

        capacity_ = old_cap * 2;
        keys_.resize(capacity_);
        values_.resize(capacity_);
        occupied_.assign(capacity_, 0);
        size_ = 0;

        size_t mask = capacity_ - 1;
        for (size_t i = 0; i < old_cap; ++i) {
            if (old_occupied[i]) {
                size_t idx = hasher_(old_keys[i]) & mask;
                while (occupied_[idx]) {
                    idx = (idx + 1) & mask;
                }
                occupied_[idx] = 1;
                keys_[idx] = std::move(old_keys[i]);
                values_[idx] = std::move(old_values[i]);
                size_++;
            }
        }
    }

    size_t capacity_ = 0;
    size_t size_ = 0;
    std::vector<Key> keys_;
    std::vector<Value> values_;
    std::vector<uint8_t> occupied_;
    Hasher hasher_;
};

    // Memoization table for full polynomials P_{x, w}(q)
    FastFlatMap<PermPair, Polynomial, FastPermPairHasher> cache; 
    
    // Memoization table for leading mu coefficients mu(x, w)
    FastFlatMap<PermPair, long long, FastPermPairHasher> mu_cache;

    // --- Core Recurrences ---
    Polynomial get_kl_recursive(Permutation x, Permutation w, bool already_bruhat_checked);
    long long get_mu_recursive(Permutation x, Permutation w, bool already_bruhat_checked);

    // --- Correction Terms for P_{x, w}(q) ---
    void apply_correction_term(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void find_extremal_pairs_and_correct(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void find_covering_relations_and_correct(Polynomial& poly, const Permutation& x, const Permutation& ws, int s_index);
    void perform_correction(Polynomial& poly, const Permutation& x, const Permutation& ws, const Permutation& z);

    // --- Correction Terms for mu(x, w) ---
    void mu_apply_correction_term(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_find_extremal_pairs_and_correct(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_find_covering_relations_and_correct(long long& current_mu, const Permutation& x, const Permutation& ws, int s_index);
    void mu_perform_correction(long long& current_mu, const Permutation& x, const Permutation& ws, const Permutation& z);
};
