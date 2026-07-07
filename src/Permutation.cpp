#include "Permutation.hpp"
#include <numeric>
#include <algorithm>
#include <stdexcept>

#include <bit>

Permutation::Permutation(const ValueType* data_ptr, int n_val) {
    if (n_val > MAX_SIZE) throw std::length_error("Permutation size exceeds maximum allowed size");
    n = n_val;
    // Copy directly without casting since the underlying type matches
    std::copy(data_ptr, data_ptr + n, data.begin());
}

Permutation::Permutation(const int* data_ptr, int n_val) {
    if (n_val > MAX_SIZE) throw std::length_error("Permutation size exceeds maximum allowed size");
    n = n_val;
    for(int i = 0; i < n; ++i) {
        data[i] = static_cast<ValueType>(data_ptr[i]);
    }
}

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

int Permutation::length() const {
    // Return cached length if it has already been computed
    if (length_ != -1) {
        return length_;
    }

    int inversions = 0;
    if (n <= 32) {
        // O(n) bitmask algorithm for N <= 32
        uint32_t seen = 0;
        for (int i = 0; i < n; ++i) {
            int v = data[i];
            if (v < 31) {
                inversions += std::popcount(seen >> (v + 1));
            }
            seen |= (1u << v);
        }
    } else {
        // Fallback to standard scalar calculation if N > 32
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (data[i] > data[j]) {
                    inversions++;
                }
            }
        }
    }
    
    length_ = inversions;
    return length_;
}

PackedInt Permutation::pack() const {
#ifdef USE_128BIT_PACKING
    if (n <= 12) {
        uint64_t low = 0;
        for (int i = 0; i < n; ++i) {
            low |= ((uint64_t)(data[i] & PACKED_MASK) << (BITS_PER_ELEMENT * i));
        }
        return low;
    }
#endif
    PackedInt packed = 0;
    for (int i = 0; i < n; ++i) {
        packed |= ((PackedInt)(data[i] & PACKED_MASK) << (BITS_PER_ELEMENT * i));
    }
    return packed;
}

std::vector<int> Permutation::unpack_to_vector(PackedInt packed, int n) {
    std::vector<int> d(n);
    for (int i = 0; i < n; ++i) {
        d[i] = (packed >> (BITS_PER_ELEMENT * i)) & PACKED_MASK;
    }
    return d;
}

Permutation Permutation::from_packed(PackedInt packed, int n) {
    Permutation p;
    p.set_size(n);
    for (int i = 0; i < n; ++i) {
        p.data[i] = static_cast<ValueType>((packed >> (BITS_PER_ELEMENT * i)) & PACKED_MASK);
    }
    return p;
}

Permutation Permutation::inverse() const {
    Permutation new_perm;
    new_perm.n = n;
    
    for (int i = 0; i < n; ++i) {
        new_perm.data[static_cast<int>(data[i])] = static_cast<ValueType>(i);
    }
    
    // The length of a permutation is always equal to the length of its inverse
    if (this->length_ != -1) {
        new_perm.length_ = this->length_;
    }
    return new_perm;
}

Permutation Permutation::multiply_simple_transposition(int i) const {
    if (i < 0 || i >= n - 1) {
        throw std::out_of_range("Index for simple transposition is out of range.");
    }
    // Perform an efficient copy of the fixed-size array
    Permutation new_perm = *this;
    
    std::swap(new_perm.data[i], new_perm.data[i + 1]);
    
    // By default, invalidate the length cache for the new permutation
    new_perm.length_ = -1;

    // If the length of the current permutation is known, we can calculate
    // the new length in O(1) time based on whether we swapped an ascent or descent.
    if (this->length_ != -1) {
        if (data[i] < data[i+1]) {
            new_perm.length_ = this->length_ + 1;
        } else {
            new_perm.length_ = this->length_ - 1;
        }
    }
    return new_perm;
}

bool Permutation::operator==(const Permutation& rhs) const {
    if (n != rhs.n) return false;
    
    // Compare elements sequentially to avoid issues with uninitialized padding
    for(int i = 0; i < n; ++i) {
        if (data[i] != rhs.data[i]) return false;
    }
    return true;
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
