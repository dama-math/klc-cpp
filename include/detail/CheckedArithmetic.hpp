#pragma once

#include <limits>
#include <stdexcept>

namespace klc::detail {

// Overflow must be reported rather than silently producing an incorrect polynomial.
inline long long checked_add(long long a, long long b) {
    long long result;
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_add_overflow(a, b, &result))
        throw std::overflow_error("Polynomial coefficient overflow");
#else
    if ((b > 0 && a > std::numeric_limits<long long>::max() - b) ||
        (b < 0 && a < std::numeric_limits<long long>::min() - b))
        throw std::overflow_error("Polynomial coefficient overflow");
    result = a + b;
#endif
    return result;
}

inline long long checked_subtract(long long a, long long b) {
    long long result;
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_sub_overflow(a, b, &result))
        throw std::overflow_error("Polynomial coefficient overflow");
#else
    if ((b > 0 && a < std::numeric_limits<long long>::min() + b) ||
        (b < 0 && a > std::numeric_limits<long long>::max() + b))
        throw std::overflow_error("Polynomial coefficient overflow");
    result = a - b;
#endif
    return result;
}

inline long long checked_multiply(long long a, long long b) {
    long long result;
#if defined(__GNUC__) || defined(__clang__)
    if (__builtin_mul_overflow(a, b, &result))
        throw std::overflow_error("Polynomial coefficient overflow");
#else
    constexpr auto lo = std::numeric_limits<long long>::min();
    constexpr auto hi = std::numeric_limits<long long>::max();
    if (a != 0 && b != 0 &&
        ((a > 0 && b > 0 && a > hi / b) ||
         (a > 0 && b < 0 && b < lo / a) ||
         (a < 0 && b > 0 && a < lo / b) ||
         (a < 0 && b < 0 && a < hi / b)))
        throw std::overflow_error("Polynomial coefficient overflow");
    result = a * b;
#endif
    return result;
}

} // namespace klc::detail
