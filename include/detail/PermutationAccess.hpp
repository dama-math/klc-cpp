#pragma once

#include "Permutation.hpp"

namespace klc::detail {

// Only reduction and candidate enumeration may temporarily mutate raw storage.
// Call commit after a mutation; public constructors always validate input.
struct PermutationAccess {
    static auto* buffer(Permutation& p) { return p.data.data(); }
    static int cached_length(const Permutation& p) { return p.length_; }

    static void commit(Permutation& p, int size, int length = -1) {
        assert(size >= 0 && size <= Permutation::MAX_SIZE);
        p.n = size;
        p.length_ = length;
#ifndef NDEBUG
        const Permutation checked(p.data.data(), size);
        assert(length == -1 || length == checked.length());
#endif
    }

    static Permutation candidate(const Permutation::ValueType* data, int size, int length) {
        Permutation p;
        std::copy_n(data, size, p.data.begin());
        commit(p, size, length);
        return p;
    }

    static Permutation covering_swap(const Permutation& p, int i, int j) {
        Permutation result = p;
        std::swap(result.data[i], result.data[j]);
        commit(result, p.size(), p.length() - 1);
        return result;
    }
};

} // namespace klc::detail
