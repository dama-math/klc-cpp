#pragma once

#include "detail/RankOperations.hpp"
#include "detail/PermutationAccess.hpp"
#include <bit>

namespace klc::detail {

// One enumeration algorithm, two rank backends. Dispatch happens once per search.
template<class Rank, class Callback>
void enumerate_extremal_pairs(const Permutation& x, const Permutation& ws,
                             int s_index, Callback&& callback) {
    const int n = x.size();
    if (n == 0) return;
    if (n == 1) { callback(x); return; }
    const auto lower = Rank::matrix(x);
    const auto upper = Rank::matrix(ws);
    const uint32_t right_descents = descent_mask(ws) | (uint32_t{1} << s_index);
    const uint32_t left_descents = descent_mask(ws.inverse());
    std::array<Permutation::ValueType, Permutation::MAX_SIZE> path{};

    auto dfs = [&](auto&& self, int depth, uint32_t used, int length,
                   typename Rank::Row column) -> bool {
        int max_value = n - 1;
        if (depth > 0 && (right_descents & (uint32_t{1} << (depth - 1)))) {
            if (path[depth - 1] == 0) return true;
            max_value = path[depth - 1] - 1;
        }
        const uint32_t max_mask = (uint32_t{1} << (max_value + 1)) - 1;
        const uint32_t rejected = left_descents & (~used >> 1);
        uint32_t available = ~used & max_mask & ~rejected;
        while (available) {
            const int value = std::countr_zero(available);
            available &= available - 1;
            const auto next = Rank::extend(column, value);
            if (!Rank::within_bounds(next, lower[depth], upper[depth])) continue;
            path[depth] = static_cast<Permutation::ValueType>(value);
            const int next_length = length + std::popcount(used >> (value + 1));
            if (depth == n - 1) {
                auto z = PermutationAccess::candidate(path.data(), n, next_length);
                if (!callback(z)) return false;
            } else if (!self(self, depth + 1, used | (uint32_t{1} << value),
                             next_length, next)) return false;
        }
        return true;
    };
    dfs(dfs, 0, 0, 0, Rank::zero());
}

template<class Callback>
void search_extremal_pairs(const Permutation& x, const Permutation& ws,
                          int s_index, Callback&& callback) {
    assert(x.size() == ws.size());
    assert(x.size() <= 1 || (s_index >= 0 && s_index + 1 < x.size()));
#if KLC_HAS_SIMD
    if (x.size() <= SimdRank::width) {
        enumerate_extremal_pairs<SimdRank>(x, ws, s_index, callback);
        return;
    }
#endif
    enumerate_extremal_pairs<ScalarRank>(x, ws, s_index, callback);
}

} // namespace klc::detail
