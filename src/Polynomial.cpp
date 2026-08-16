#include "Polynomial.hpp"
#include <algorithm>
#include <stdexcept>

Polynomial::Polynomial(std::initializer_list<long long> coeffs) {
    deg = static_cast<int>(coeffs.size()) - 1;
    if (deg >= MAX_CAPACITY) {
        throw std::length_error("Polynomial degree exceeds MAX_CAPACITY");
    }
    data.fill(0);
    size_t i = 0;
    for (long long c : coeffs) {
        data[i++] = c;
    }
    trim();
}

Polynomial::Polynomial(const std::vector<long long>& coeffs) {
    deg = static_cast<int>(coeffs.size()) - 1;
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

std::string Polynomial::to_string_simple() const {
    if (is_zero()) return "0";
    
    std::string s;
    for (int i = 0; i <= deg; ++i) {
        if (i > 0) s += ',';
        s += std::to_string(data[i]);
    }
    return s;
}
