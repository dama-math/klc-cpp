#pragma once

#include <vector>
#include <string>
#include <numeric>
#include <array>
#include <cstdint>

/**
 * @brief Represents a polynomial.
 * 
 * Highly efficient fixed-size polynomial representation.
 * Avoids any heap allocation or variant overhead.
 */
class Polynomial {
public:
    static constexpr int MAX_CAPACITY = 32;

    /** @brief Default constructor creating a zero polynomial. */
    Polynomial() : deg(-1) {}

    /** @brief Constructs a polynomial from a vector of coefficients. */
    Polynomial(const std::vector<long long>& coeffs);

    /** @return The degree of the polynomial. */
    int degree() const { return deg; }

    /** @return Coefficient of the i-th power of q. */
    long long get_coeff(int i) const {
        if (i < 0 || i > deg) return 0;
        return data[i];
    }

    /** @brief Equality operator. */
    bool operator==(const Polynomial& rhs) const;

    /** @brief Addition assignment operator. */
    Polynomial& operator+=(const Polynomial& rhs);

    /** @brief Subtraction assignment operator. */
    Polynomial& operator-=(const Polynomial& rhs);

    /** @brief Addition operator. */
    Polynomial operator+(const Polynomial& rhs) const;

    /** @brief Subtraction operator. */
    Polynomial operator-(const Polynomial& rhs) const;

    /** @brief Multiplication assignment operator. */
    Polynomial& operator*=(const Polynomial& rhs);

    /** @brief Multiplication operator. */
    Polynomial operator*(const Polynomial& rhs) const;

    /**
     * @brief Shifts the polynomial by multiplying by q^amount.
     * @param amount The power of q to multiply by.
     * @return Shifted polynomial.
     */
    Polynomial shift(int amount) const;

    /**
     * @brief Multiplies all coefficients by a scalar value.
     * @param scalar The scalar value.
     * @return Scaled polynomial.
     */
    Polynomial multiply_by_scalar(long long scalar) const;

    /** @brief Directly adds rhs * q^shift to *this in-place without creating temporaries. */
    void add_shifted(const Polynomial& rhs, int shift, long long scalar = 1) {
        if (rhs.deg < 0 || scalar == 0) return;
        int needed_deg = rhs.deg + shift;
        if (deg < needed_deg) {
            std::fill(data.begin() + deg + 1, data.begin() + needed_deg + 1, 0);
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

    /** @brief Directly subtracts rhs * q^shift to *this in-place without creating temporaries. */
    void subtract_shifted(const Polynomial& rhs, int shift, long long scalar = 1) {
        if (rhs.deg < 0 || scalar == 0) return;
        int needed_deg = rhs.deg + shift;
        if (deg < needed_deg) {
            std::fill(data.begin() + deg + 1, data.begin() + needed_deg + 1, 0);
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

    /**
     * @brief Evaluates the polynomial at a given value of q.
     * @param q The value to substitute for q.
     * @return Evaluated result.
     */
    long long substitute(long long q) const;

    /** @return String representation of the polynomial. */
    std::string to_string() const;

    /** @return Simplified string representation (compact). */
    std::string to_string_simple() const;

    /** @return True if the polynomial is identically zero. */
    bool is_zero() const { return deg < 0; }

    const long long* data_ptr() const { return data.data(); }
    long long* data_ptr_mut() { return data.data(); }

private:
    int deg = -1; // -1 for zero polynomial
    std::array<long long, MAX_CAPACITY> data{};

    void trim() {
        while (deg >= 0 && data[deg] == 0) {
            deg--;
        }
    }
};


