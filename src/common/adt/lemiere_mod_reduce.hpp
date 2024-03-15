//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdint>

namespace mod_reduce {

// Utility function to compute (x * y) >> 64, or "multiply high".
// On x86-64, this is a single instruction, but not all platforms
// support the __uint128_t type, so we provide a generic
// implementation as well.
inline constexpr uint64_t multiply_high_u64(uint64_t x, uint64_t y) {
#if defined(__SIZEOF_INT128__)
    return (uint64_t)(((__uint128_t)x * (__uint128_t)y) >> 64);
#else
    // For platforms without int128 support, do it the long way.
    uint64_t x_lo = x & 0xffffffff;
    uint64_t x_hi = x >> 32;
    uint64_t buckets_lo = y & 0xffffffff;
    uint64_t buckets_hi = y >> 32;
    uint64_t prod_hi = x_hi * buckets_hi;
    uint64_t prod_lo = x_lo * buckets_lo;
    uint64_t prod_mid1 = x_hi * buckets_lo;
    uint64_t prod_mid2 = x_lo * buckets_hi;
    uint64_t carry =
            ((prod_mid1 & 0xffffffff) + (prod_mid2 & 0xffffffff) + (prod_lo >> 32)) >>
            32;
    return prod_hi + (prod_mid1 >> 32) + (prod_mid2 >> 32) + carry;
#endif
}

}

