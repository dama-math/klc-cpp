#pragma once

#include <array>
#include <initializer_list>
#include <string>
#include <vector>
#include <cstdint>
#include <algorithm>
#include <cstring>

/**
 * @brief Represents a single-variable polynomial in q with integer coefficients.
 * 
 * Uses a fixed-size array buffer to eliminate heap allocations.
 * Degree `deg = -1` represents the zero polynomial.
 */
class Polynomial {
public:
    static constexpr int MAX_CAPACITY = 64;

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
        if (rhs.deg < 0 || scalar == 0) return;
        int needed_deg = rhs.deg + shift;
        if (deg < needed_deg) {
            for (int i = deg + 1; i <= needed_deg; ++i) data[i] = 0;
            deg = needed_deg;
        }
        if (scalar == 1) {
            for (int i = 0; i <= rhs.deg; ++i) {
                data[i + shift] += rhs.data[i];
            }
        } else {
            for (int i = 0; i <= rhs.deg; ++i) {
                data[i + shift] += rhs.data[i] * scalar;
            }
        }
        trim();
    }

    /**
     * @brief In-place subtraction: *this -= rhs * (scalar * q^shift).
     * Avoids creating temporary objects during recursive correction subtraction.
     */
    void subtract_shifted(const Polynomial& rhs, int shift, long long scalar = 1) {
        if (rhs.deg < 0 || scalar == 0) return;
        int needed_deg = rhs.deg + shift;
        if (deg < needed_deg) {
            for (int i = deg + 1; i <= needed_deg; ++i) data[i] = 0;
            deg = needed_deg;
        }
        if (scalar == 1) {
            for (int i = 0; i <= rhs.deg; ++i) {
                data[i + shift] -= rhs.data[i];
            }
        } else {
            for (int i = 0; i <= rhs.deg; ++i) {
                data[i + shift] -= rhs.data[i] * scalar;
            }
        }
        trim();
    }

    /** @brief Evaluates the polynomial at q (e.g. q = 1). */
    long long substitute(long long q) const;

    /** @brief Formats coefficients as comma-separated values: "a_0,a_1,a_2,...". */
    std::string to_string_simple() const;

private:
    int deg = -1;
    std::array<long long, MAX_CAPACITY> data{};

    void trim() {
        while (deg >= 0 && data[deg] == 0) {
            deg--;
        }
    }
};
