#include "Polynomial.hpp"
#include <algorithm>
#include <iostream>
#include <cstring>  // for memcpy
#include <cassert>

// Constructors
Polynomial::Polynomial() : deg(-1), storage(std::in_place_type<SmallStorage>) {
    std::get<SmallStorage>(storage).fill(0);
}

Polynomial::Polynomial(const std::vector<long long>& coeffs) {
    int n = coeffs.size();
    if (n <= SMALL_SIZE) {
        storage.emplace<SmallStorage>();
        auto& arr = std::get<SmallStorage>(storage);
        arr.fill(0);
        std::copy(coeffs.begin(), coeffs.end(), arr.begin());
        deg = n - 1;
    } else {
        storage.emplace<LargeStorage>(coeffs);
        deg = n - 1;
    }
    trim();
}

// Accessors
int Polynomial::degree() const {
    return deg;
}

long long Polynomial::get_coeff(int i) const {
    if (i < 0 || i > deg) return 0;
    return data_ptr()[i];
}

std::vector<long long> Polynomial::get_coeffs() const {
    if (deg < 0) return {0};
    const long long* ptr = data_ptr();
    return std::vector<long long>(ptr, ptr + deg + 1);
}

// Operators
bool Polynomial::operator==(const Polynomial& rhs) const {
    if (deg != rhs.deg) return false;
    if (deg < 0) return true;
    
    const long long* p1 = data_ptr();
    const long long* p2 = rhs.data_ptr();
    
    return std::equal(p1, p1 + deg + 1, p2);
}

void Polynomial::resize_storage(int new_size) {
    if (new_size <= SMALL_SIZE) {
        if (std::holds_alternative<LargeStorage>(storage)) {
            SmallStorage arr;
            arr.fill(0);
            auto& vec = std::get<LargeStorage>(storage);
            std::copy(vec.begin(), vec.begin() + new_size, arr.begin());
            storage.emplace<SmallStorage>(arr);
        } else {
            auto& arr = std::get<SmallStorage>(storage);
            if (deg < new_size - 1) {
                 std::fill(arr.begin() + deg + 1, arr.begin() + new_size, 0);
            }
        }
    } else {
        if (std::holds_alternative<SmallStorage>(storage)) {
            switch_to_large();
        }
        std::get<LargeStorage>(storage).resize(new_size, 0);
    }
}

void Polynomial::switch_to_large() {
    if (std::holds_alternative<LargeStorage>(storage)) return;
    auto& arr = std::get<SmallStorage>(storage);
    LargeStorage vec(std::max((int)deg + 1, (int)SMALL_SIZE), 0);
    std::copy(arr.begin(), arr.begin() + deg + 1, vec.begin());
    storage.emplace<LargeStorage>(std::move(vec));
}

Polynomial& Polynomial::operator+=(const Polynomial& rhs) {
    if (rhs.deg < 0) return *this;
    
    int new_deg = std::max(deg, rhs.deg);
    int current_capacity = std::holds_alternative<SmallStorage>(storage) ? SMALL_SIZE : std::get<LargeStorage>(storage).size();
    if (new_deg + 1 > current_capacity) {
        resize_storage(new_deg + 1);
    }
    
    long long* p1 = data_ptr_mut();
    const long long* p2 = rhs.data_ptr();
    
    for (int i = 0; i <= rhs.deg; ++i) {
        p1[i] += p2[i];
    }
    
    deg = new_deg;
    trim();
    return *this;
}

Polynomial& Polynomial::operator-=(const Polynomial& rhs) {
    if (rhs.deg < 0) return *this;

    int new_deg = std::max(deg, rhs.deg);
    int current_capacity = std::holds_alternative<SmallStorage>(storage) ? SMALL_SIZE : std::get<LargeStorage>(storage).size();
    if (new_deg + 1 > current_capacity) {
        resize_storage(new_deg + 1);
    }

    long long* p1 = data_ptr_mut();
    const long long* p2 = rhs.data_ptr();

    for (int i = 0; i <= rhs.deg; ++i) {
        p1[i] -= p2[i];
    }

    deg = new_deg;
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
    if (this->is_zero() || rhs.is_zero()) {
        deg = -1;
        return *this;
    }

    // Guard against self-multiplication (p *= p) which would invalidate
    // the rhs data pointer when storage is modified in-place.
    if (this == &rhs) {
        Polynomial copy(rhs);
        return *this *= copy;
    }
    
    int new_deg = deg + rhs.deg;
    int needed_size = new_deg + 1;
    
    if (needed_size <= SMALL_SIZE) {
        SmallStorage buf;
        buf.fill(0);
        
        const long long* A = data_ptr();
        const long long* B = rhs.data_ptr();
        
        for (int i = 0; i <= deg; ++i) {
            if (A[i] == 0) continue;
            for (int j = 0; j <= rhs.deg; ++j) {
                buf[i + j] += A[i] * B[j];
            }
        }
        
        storage.emplace<SmallStorage>(buf);
        deg = new_deg;
    } else {
        LargeStorage res(needed_size, 0);
        const long long* A = data_ptr();
        const long long* B = rhs.data_ptr();
        
        for (int i = 0; i <= deg; ++i) {
            if (A[i] == 0) continue;
            for (int j = 0; j <= rhs.deg; ++j) {
                res[i + j] += A[i] * B[j];
            }
        }
        
        storage.emplace<LargeStorage>(std::move(res));
        deg = new_deg;
    }
    
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
    if (amount == 0) return *this;
    if (is_zero()) return *this;

    Polynomial res;
    int new_deg = deg + amount;
    int needed_size = new_deg + 1;
    
    res.resize_storage(needed_size);
    long long* dst = res.data_ptr_mut();
    const long long* src = data_ptr();
    
    std::copy(src, src + deg + 1, dst + amount);
    
    res.deg = new_deg;
    return res;
}

Polynomial Polynomial::multiply_by_scalar(long long scalar) const {
    if (scalar == 0) return Polynomial(); 
    if (scalar == 1) return *this;
    
    Polynomial res(*this);
    long long* ptr = res.data_ptr_mut();
    for(int i = 0; i <= deg; ++i) {
        ptr[i] *= scalar;
    }
    return res;
}

long long Polynomial::substitute(long long q) const {
    if (deg < 0) return 0;
    long long result = 0;
    long long q_power = 1;
    const long long* ptr = data_ptr();
    for (int i = 0; i <= deg; ++i) {
        result += ptr[i] * q_power;
        q_power *= q;
    }
    return result;
}

std::string Polynomial::to_string() const {
    if (is_zero()) return "0";

    std::string s = "";
    const long long* coeffs = data_ptr();
    
    for (int i = deg; i >= 0; --i) {
        if (coeffs[i] != 0) {
            if (i < deg) {
                s += (coeffs[i] > 0) ? " + " : " - ";
            } else if (coeffs[i] < 0) {
                 s += "-";
            }
            
            long long abs_coeff = std::abs(coeffs[i]);
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
    const long long* coeffs = data_ptr();
    for (int i = 0; i <= deg; ++i) {
        if (i > 0) s += ",";
        s += std::to_string(coeffs[i]);
    }
    return s;
}

bool Polynomial::is_zero() const {
    return deg < 0;
}

std::size_t Polynomial::hash() const {
    std::size_t seed = 0;
    int current_deg = degree();
    seed ^= std::hash<int>{}(current_deg) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    
    if (current_deg >= 0) {
        const long long* ptr = data_ptr();
        for (int i = 0; i <= current_deg; ++i) {
            seed ^= std::hash<long long>{}(ptr[i]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    }
    return seed;
}

void Polynomial::trim() {
    long long* ptr = data_ptr_mut();
    while (deg >= 0 && ptr[deg] == 0) {
        deg--;
    }
}
