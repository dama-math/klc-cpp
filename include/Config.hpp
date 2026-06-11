#pragma once

#include <cstdint>

// Define USE_128BIT_PACKING via compiler flags (-DUSE_128BIT_PACKING) for n > 16.
// By default, if not defined, it will use 64-bit packing for max performance up to n=16.

#ifdef USE_128BIT_PACKING
    // 128-bit mode: Safe for N <= 24 (using 5 bits per element)
    using PackedInt = __uint128_t;
    constexpr int BITS_PER_ELEMENT = 5;
    constexpr int MAX_PACKED_SIZE = 24;
    constexpr PackedInt PACKED_MASK = 0x1F;
#else
    // 64-bit mode: Max performance for N <= 16 (using 4 bits per element)
    using PackedInt = uint64_t;
    constexpr int BITS_PER_ELEMENT = 4;
    constexpr int MAX_PACKED_SIZE = 16;
    constexpr PackedInt PACKED_MASK = 0xF;
#endif
