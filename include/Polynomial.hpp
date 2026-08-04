#pragma once

#include <vector>
#include <string>
#include <numeric>
#include <array>
#include <variant>

/**
 * @brief Represents a polynomial.
 * 
 * Highly efficient polynomial representation using std::variant for Small Buffer Optimization.
 * Small polynomials (degree < 8) are stored directly in the object.
 * Larger polynomials fall back to std::vector.
 */
class Polynomial {
public:
    /** @brief Default constructor creating a zero polynomial. */
    Polynomial();

    /** @brief Constructs a polynomial from a vector of coefficients. */
    Polynomial(const std::vector<long long>& coeffs);

    /** @return The degree of the polynomial. */
    int degree() const;

    /** @return Coefficient of the i-th power of q. */
    long long get_coeff(int i) const;


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
    bool is_zero() const;

private:
    static constexpr int SMALL_SIZE = 8;
    int deg = -1; // -1 for zero polynomial
    
    using SmallStorage = std::array<long long, SMALL_SIZE>;
    using LargeStorage = std::vector<long long>;
    std::variant<SmallStorage, LargeStorage> storage;

    void trim();
    void switch_to_large(); // Helper to promote to large storage
    
    // Helpers for internal access
    long long const* data_ptr() const {
        if (const auto* p = std::get_if<SmallStorage>(&storage)) return p->data();
        return std::get_if<LargeStorage>(&storage)->data();
    }
    
    long long* data_ptr_mut() {
        if (auto* p = std::get_if<SmallStorage>(&storage)) return p->data();
        return std::get_if<LargeStorage>(&storage)->data();
    }
    
    int data_size() const {
        return deg + 1;
    }
    
    void resize_storage(int new_size);
};

