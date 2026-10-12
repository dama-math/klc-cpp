#pragma once

#include <array>
#include <initializer_list>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <stdexcept>
#include "Config.hpp"
#include "detail/CheckedArithmetic.hpp"

/**
 * @brief Represents a single-variable polynomial in q with integer coefficients.
 * 
 * Uses a fixed-size array buffer to eliminate heap allocations.
 * Degree `deg = -1` represents the zero polynomial.
 */
class Polynomial {
public:
    // Retain the compact 64-bit-mode layout; accommodate S_24 in large mode.
    static constexpr int MAX_CAPACITY = MAX_PACKED_SIZE <= 16 ? 64 :
        MAX_PACKED_SIZE * (MAX_PACKED_SIZE - 1) / 4 + 1;

    /** @brief Default constructor creating a zero polynomial (degree -1). */
    Polynomial() : deg(-1) {}

    /** @brief Constructs a polynomial from a list of coefficients [a_0, a_1, ...]. */
    Polynomial(std::initializer_list<long long> coeffs);

    /** @brief Constructs a polynomial from a vector of coefficients. */
    Polynomial(const std::vector<long long>& coeffs);

    int degree() const { return deg; }
    bool is_zero() const { return deg < 0; }

    /** @brief Returns coefficient of q^i (0 if out of range). */
    long long get_coeff(int i) const {
        if (i < 0 || i > deg) return 0;
        return data[i];
    }

    bool operator==(const Polynomial& rhs) const;

    /** @brief In-place polynomial multiplication (*this *= rhs). */
    Polynomial& operator*=(const Polynomial& rhs);

    /**
     * @brief In-place addition: *this += rhs * (scalar * q^shift).
     * Avoids creating temporary objects during recursive sum accumulation.
     */
    void add_shifted(const Polynomial& rhs, int shift, long long scalar = 1) {
        accumulate_shifted<false>(rhs, shift, scalar);
    }

    /**
     * @brief In-place subtraction: *this -= rhs * (scalar * q^shift).
     * Avoids creating temporary objects during recursive correction subtraction.
     */
    void subtract_shifted(const Polynomial& rhs, int shift, long long scalar = 1) {
        accumulate_shifted<true>(rhs, shift, scalar);
    }

    /** @brief Evaluates the polynomial at q (e.g. q = 1). */
    long long substitute(long long q) const;

    /** @brief Formats coefficients as comma-separated values: "a_0,a_1,a_2,...". */
    std::string to_string_simple() const;

private:
    int deg = -1;
    std::array<long long, MAX_CAPACITY> data{};

    template<bool Subtract>
    void accumulate_shifted(const Polynomial& rhs, int shift, long long scalar) {
        if (shift < 0) throw std::invalid_argument("Polynomial shift must be nonnegative");
        if (rhs.deg < 0 || scalar == 0) return;
        if (shift >= MAX_CAPACITY || rhs.deg >= MAX_CAPACITY - shift)
            throw std::length_error("Polynomial degree exceeds MAX_CAPACITY");
        // Descending traversal preserves rhs when adding a shifted polynomial to itself.
        const int rhs_degree = rhs.deg;
        const int needed_degree = rhs_degree + shift;
        for (int i = deg + 1; i <= needed_degree; ++i) data[i] = 0;
        for (int i = rhs_degree; i >= 0; --i) {
            const auto term = scalar == 1 ? rhs.data[i] :
                klc::detail::checked_multiply(rhs.data[i], scalar);
            if constexpr (Subtract)
                data[i + shift] = klc::detail::checked_subtract(data[i + shift], term);
            else
                data[i + shift] = klc::detail::checked_add(data[i + shift], term);
        }
        deg = std::max(deg, needed_degree);
        trim();
    }

    void trim() {
        while (deg >= 0 && data[deg] == 0) {
            deg--;
        }
    }
};
