#pragma once

#include "Permutation.hpp"
#include <array>
#include <cstdint>

// The optimized backend requires SSE4.1; other targets use the scalar backend.
#if defined(__SSE4_1__) && !defined(KLC_DISABLE_SIMD)
#define KLC_HAS_SIMD 1
#include <immintrin.h>
#else
#define KLC_HAS_SIMD 0
#endif

namespace klc::detail {

inline uint32_t descent_mask(const Permutation& p) {
    uint32_t mask = 0;
    for (int i = 0; i + 1 < p.size(); ++i)
        if (p[i] > p[i + 1]) mask |= uint32_t{1} << i;
    return mask;
}

struct ScalarRank {
    using Row = std::array<int8_t, Permutation::MAX_SIZE>;
    using Matrix = std::array<Row, Permutation::MAX_SIZE>;
    static Row zero() { return {}; }
    static Row extend(Row row, int value) {
        for (int k = 0; k <= value; ++k) ++row[k];
        return row;
    }
    static bool within_bounds(const Row& row, const Row& lower, const Row& upper) {
        for (int k = 0; k < Permutation::MAX_SIZE; ++k)
            if (row[k] < lower[k] || row[k] > upper[k]) return false;
        return true;
    }
    static Matrix matrix(const Permutation& p) {
        Matrix result{};
        Row row{};
        for (int i = 0; i < p.size(); ++i) {
            row = extend(row, p[i]);
            result[i] = row;
        }
        return result;
    }
};

#if KLC_HAS_SIMD
struct SimdRank {
    static constexpr int width = 16;
    struct alignas(16) Row { __m128i value; };
    using Matrix = std::array<Row, width>;

    struct alignas(16) Indicators {
        int8_t data[width][width]{};
        constexpr Indicators() {
            for (int i = 0; i < width; ++i)
                for (int k = 0; k <= i; ++k) data[i][k] = 1;
        }
    };
    static Row zero() { return {_mm_setzero_si128()}; }
    static Row extend(Row row, int value) {
        assert(value >= 0 && value < width);
        static constexpr Indicators indicators;
        auto mask = _mm_load_si128(reinterpret_cast<const __m128i*>(indicators.data[value]));
        return {_mm_add_epi8(row.value, mask)};
    }
    static bool within_bounds(Row row, Row lower, Row upper) {
        auto bad = _mm_or_si128(_mm_cmpgt_epi8(row.value, upper.value),
                               _mm_cmpgt_epi8(lower.value, row.value));
        return _mm_testz_si128(bad, bad);
    }
    static Matrix matrix(const Permutation& p) {
        assert(p.size() <= width);
        Matrix result{};
        auto row = zero();
        for (int i = 0; i < p.size(); ++i) {
            row = extend(row, p[i]);
            result[i] = row;
        }
        return result;
    }
};
#endif

inline bool scalar_bruhat_compare(const Permutation& x, const Permutation& w) {
    std::array<int8_t, Permutation::MAX_SIZE> diff{};
    for (int i = 0; i < x.size(); ++i) {
        const int a = x[i], b = w[i];
        if (a < b) {
            for (int k = a; k < b; ++k) ++diff[k];
        } else {
            for (int k = b; k < a; ++k)
                if (--diff[k] < 0) return false;
        }
    }
    return true;
}

#if KLC_HAS_SIMD
inline bool simd_bruhat_compare(const Permutation& x, const Permutation& w) {
    struct alignas(16) Intervals {
        int8_t data[16][16][16]{};
        constexpr Intervals() {
            for (int a = 0; a < 16; ++a)
                for (int b = a; b < 16; ++b)
                    for (int k = a; k < b; ++k) data[a][b][k] = 1;
        }
    };
    static constexpr Intervals intervals;
    auto diff = _mm_setzero_si128();
    for (int i = 0; i < x.size(); ++i) {
        const int a = x[i], b = w[i];
        if (a < b) {
            auto mask = _mm_load_si128(reinterpret_cast<const __m128i*>(intervals.data[a][b]));
            diff = _mm_add_epi8(diff, mask);
        } else if (a > b) {
            auto mask = _mm_load_si128(reinterpret_cast<const __m128i*>(intervals.data[b][a]));
            diff = _mm_sub_epi8(diff, mask);
            if (_mm_movemask_epi8(diff)) return false;
        }
    }
    return true;
}
#endif

} // namespace klc::detail
