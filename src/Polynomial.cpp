#include "Polynomial.hpp"
#include <algorithm>
#include <iostream>
#include <cstring>
#include <cassert>

Polynomial::Polynomial(const std::vector<long long>& coeffs) {
    deg = (int)coeffs.size() - 1;
    if (deg >= MAX_CAPACITY) {
        throw std::length_error("Polynomial degree exceeds MAX_CAPACITY");
    }
    data.fill(0);
    for (size_t i = 0; i < coeffs.size(); ++i) {
        data[i] = coeffs[i];
    }
    trim();
}

bool Polynomial::operator==(const Polynomial& rhs) const {
    if (deg != rhs.deg) return false;
    if (deg < 0) return true;
    return std::equal(data.begin(), data.begin() + deg + 1, rhs.data.begin());
}

Polynomial& Polynomial::operator+=(const Polynomial& rhs) {
    if (rhs.deg < 0) return *this;
    if (deg < rhs.deg) {
        std::fill(data.begin() + deg + 1, data.begin() + rhs.deg + 1, 0);
        deg = rhs.deg;
    }
    for (int i = 0; i <= rhs.deg; ++i) {
        data[i] += rhs.data[i];
    }
    trim();
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& rhs) {
    if (rhs.deg < 0) return *this;
    if (deg < rhs.deg) {
        std::fill(data.begin() + deg + 1, data.begin() + rhs.deg + 1, 0);
        deg = rhs.deg;
    }
    for (int i = 0; i <= rhs.deg; ++i) {
        data[i] -= rhs.data[i];
    }
    trim();
    return *this;
}

Polynomial Polynomial::operator+(const Polynomial& rhs) const {
    Polynomial result(*this);
    result += rhs;
    return result;
}

Polynomial Polynomial::operator-(const Polynomial& rhs) const {
    Polynomial result(*this);
    result -= rhs;
    return result;
}

Polynomial& Polynomial::operator*=(const Polynomial& rhs) {
    if (deg < 0 || rhs.deg < 0) {
        deg = -1;
        return *this;
    }
    int new_deg = deg + rhs.deg;
    if (new_deg >= MAX_CAPACITY) {
        throw std::length_error("Product degree exceeds MAX_CAPACITY");
    }
    std::array<long long, MAX_CAPACITY> buf{};
    for (int i = 0; i <= deg; ++i) {
        long long c = data[i];
        if (c == 0) continue;
        for (int j = 0; j <= rhs.deg; ++j) {
            buf[i + j] += c * rhs.data[j];
        }
    }
    data = buf;
    deg = new_deg;
    trim();
    return *this;
}

Polynomial Polynomial::operator*(const Polynomial& rhs) const {
    Polynomial result(*this);
    result *= rhs;
    return result;
}

Polynomial Polynomial::shift(int amount) const {
    assert(amount >= 0 && "shift amount must be non-negative");
    if (amount == 0 || deg < 0) return *this;

    Polynomial res;
    res.deg = deg + amount;
    if (res.deg >= MAX_CAPACITY) {
        throw std::length_error("Shifted degree exceeds MAX_CAPACITY");
    }
    std::fill(res.data.begin(), res.data.begin() + amount, 0);
    std::memcpy(res.data.data() + amount, data.data(), (deg + 1) * sizeof(long long));
    return res;
}

Polynomial Polynomial::multiply_by_scalar(long long scalar) const {
    if (scalar == 0 || deg < 0) return Polynomial(); 
    if (scalar == 1) return *this;
    
    Polynomial res(*this);
    for (int i = 0; i <= deg; ++i) {
        res.data[i] *= scalar;
    }
    return res;
}

long long Polynomial::substitute(long long q) const {
    if (deg < 0) return 0;
    long long result = 0;
    long long q_power = 1;
    for (int i = 0; i <= deg; ++i) {
        result += data[i] * q_power;
        q_power *= q;
    }
    return result;
}

std::string Polynomial::to_string() const {
    if (is_zero()) return "0";

    std::string s = "";
    for (int i = deg; i >= 0; --i) {
        if (data[i] != 0) {
            if (i < deg) {
                s += (data[i] > 0) ? " + " : " - ";
            } else if (data[i] < 0) {
                 s += "-";
            }
            
            long long abs_coeff = std::abs(data[i]);
            if (abs_coeff != 1 || i == 0) {
                s += std::to_string(abs_coeff);
            }

            if (i > 0) {
                s += "q";
                if (i > 1) {
                    s += "^" + std::to_string(i);
                }
            }
        }
    }
    return s;
}

std::string Polynomial::to_string_simple() const {
    if (is_zero()) return "0";
    
    std::string s = "";
    for (int i = 0; i <= deg; ++i) {
        if (i > 0) s += ",";
        s += std::to_string(data[i]);
    }
    return s;
}

