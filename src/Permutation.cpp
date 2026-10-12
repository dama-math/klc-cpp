#include "Permutation.hpp"
#include <stdexcept>
#include <vector>

namespace {
void validate_size(int size) {
    if (size < 0 || size > Permutation::MAX_SIZE)
        throw std::length_error("Permutation size exceeds maximum allowed size (" +
                                std::to_string(Permutation::MAX_SIZE) + ")");
}

void validate_value(int value, int size, uint32_t& seen) {
    if (value < 0 || value >= size || (seen & (uint32_t{1} << value)))
        throw std::invalid_argument("Invalid permutation: elements must be a permutation of 0 to n-1");
    seen |= uint32_t{1} << value;
}
} // namespace

Permutation::Permutation(const ValueType* values, int size) {
    validate_size(size);
    if (size != 0 && values == nullptr)
        throw std::invalid_argument("Null permutation buffer");
    n = size;
    uint32_t seen = 0;
    for (int i = 0; i < n; ++i) {
        validate_value(values[i], n, seen);
        data[i] = values[i];
    }
}

Permutation Permutation::identity(int size) {
    validate_size(size);
    Permutation p;
    p.n = size;
    p.length_ = 0;
    for (int i = 0; i < size; ++i) p.data[i] = static_cast<ValueType>(i);
    return p;
}

Permutation::Permutation(const std::vector<int>& data_vec) {
    if (data_vec.size() > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = static_cast<int>(data_vec.size());
    uint32_t seen = 0;
    for (int i = 0; i < n; ++i) {
        validate_value(data_vec[i], n, seen);
        data[i] = static_cast<ValueType>(data_vec[i]);
    }
}

Permutation::Permutation(const std::string& one_line_notation) {
    parse_one_line_notation(one_line_notation);
}

void Permutation::parse_one_line_notation(const std::string& notation) {
    if (notation.size() > MAX_SIZE) {
        throw std::length_error("Permutation size exceeds maximum allowed size (" + std::to_string(MAX_SIZE) + ")");
    }
    n = static_cast<int>(notation.size());

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

        validate_value(val, n, seen_mask);
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
