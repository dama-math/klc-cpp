#include "Permutation.hpp"
#include <stdexcept>
#include <vector>

Permutation::Permutation(const std::vector<int>& data_vec) {
    if (data_vec.size() > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = static_cast<int>(data_vec.size());
    for (int i = 0; i < n; ++i) {
        data[i] = static_cast<ValueType>(data_vec[i]);
    }
}

Permutation::Permutation(const std::string& one_line_notation) {
    parse_one_line_notation(one_line_notation);
}

void Permutation::parse_one_line_notation(const std::string& notation) {
    int len = static_cast<int>(notation.size());
    if (len > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = len;

    uint32_t seen_mask = 0;
    for (int i = 0; i < n; ++i) {
        char c = notation[i];
        int val;
        if (c >= '0' && c <= '9') {
            val = c - '0';
        } else if (c >= 'a' && c <= 'z') {
            val = c - 'a' + 10;
        } else {
            throw std::invalid_argument(std::string("Invalid character in permutation string: ") + c);
        }

        if (val < 0 || val >= n || ((seen_mask >> val) & 1u)) {
            throw std::invalid_argument("Invalid permutation: elements must be a permutation of 0 to n-1");
        }
        seen_mask |= (1u << val);
        data[i] = static_cast<ValueType>(val);
    }
}

std::string Permutation::to_string() const {
    std::string s;
    s.reserve(n);
    for (int i = 0; i < n; ++i) {
        int val = static_cast<int>(data[i]);
        if (val < 10) {
            s += static_cast<char>('0' + val);
        } else {
            s += static_cast<char>('a' + (val - 10));
        }
    }
    return s;
}

std::vector<int> Permutation::to_vector() const {
    std::vector<int> res(n);
    for (int i = 0; i < n; ++i) {
        res[i] = static_cast<int>(data[i]);
    }
    return res;
}
