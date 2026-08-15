#include "Permutation.hpp"
#include <numeric>
#include <algorithm>
#include <stdexcept>
#include <bit>

Permutation::Permutation(const std::string& one_line_notation) {
    std::vector<int> vec = from_one_line_notation(one_line_notation);
    if (vec.size() > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = vec.size();
    for(int i = 0; i < n; ++i) {
        data[i] = static_cast<ValueType>(vec[i]);
    }
}

Permutation::Permutation(const std::vector<int>& data_vec) {
    if (data_vec.size() > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = data_vec.size();
    for(int i = 0; i < n; ++i) {
        data[i] = static_cast<ValueType>(data_vec[i]);
    }
}

std::vector<int> Permutation::to_vector() const {
    std::vector<int> res(n);
    for(int i = 0; i < n; ++i) {
        res[i] = static_cast<int>(data[i]);
    }
    return res;
}

std::string Permutation::to_string() const {
    return to_one_line_notation(to_vector());
}

std::string Permutation::to_one_line_notation(const std::vector<int>& data) {
    std::string s = "";
    for (int val : data) {
        if (val < 10) {
            s += std::to_string(val);
        } else {
            // Use letters for values 10 and above
            s += (char)('a' + (val - 10));
        }
    }
    return s;
}

std::vector<int> Permutation::from_one_line_notation(const std::string& notation) {
    std::vector<int> d;
    for (char c : notation) {
        if (c >= '0' && c <= '9') {
            d.push_back(c - '0');
        } else if (c >= 'a' && c <= 'z') {
            d.push_back(c - 'a' + 10);
        } else {
            throw std::invalid_argument(std::string("Invalid character in permutation string: ") + c);
        }
    }
    
    int n = d.size();
    if (n > MAX_SIZE) {
        throw std::length_error("Representation exceeds maximum size in from_one_line_notation");
    }

    std::vector<bool> seen(n, false);
    for (int val : d) {
        if (val < 0 || val >= n || seen[val]) {
            throw std::invalid_argument("Invalid permutation: not a valid permutation of 0 to n-1");
        }
        seen[val] = true;
    }
    return d;
}
