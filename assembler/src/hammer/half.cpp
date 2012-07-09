//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "half.hpp"
#include <stdint.h>

#include "toFloat.inl"

#define __MASK(size) ((size) == 32 ? 0xFFFFFFFFu : (1u << (size)) - 1)

uint16_t prob_half::convert(uint32_t i) {
  // Split stuff into sign, exponent and significand
  unsigned s =  (i >> 31);
  int      e = ((i >> 23) & 0x000000ff) - 127 + 31;
  unsigned m =   i & 0x007fffff;
  half_holder h;

  // Negative numbers are rounded to zeros
  if (s)
    return 0;

  if (e <= 0) {
    if (e < -11) {
      // Too small value, round to zero
      return 0;
    }

    // Ok, we know that e is between 0 and -11. Here i is normalized float. Try
    // to map to denormalized prob_half.
    m |= (1 << 23);

    unsigned shift = 12 - e;
    m = (m + __MASK(shift - 1) + ((m >> shift) & 1));
    m >>= shift;

    // During rounding everything might become normalized number. Handle this as
    // well.
    if (m >> 11) {
      h.m = 0;
      h.e = 1;
    } else {
      h.m = m & __MASK(11);
      h.e = 0;
    }
  } else {
    // i is normalized floating point which will map to normalized prob_half
    // (possible infinite).

    // Handle exponent overflow.
    if (e > 31)
      return 0xFFFF;

    // Simple case (normal number) - round the significand, m, to 11
    // bits and combine it with exponent.
    h.e = e;
    h.m = (m + __MASK(23 - 11 - 1) + ((m >> (23 - 11)) & 1)) >> (23 - 11);
    return h.bits;
  }

  return h.bits;
}
