//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _PROB_HALF_H_
#define _PROB_HALF_H_

#include <iostream>

class prob_half {
public:

  prob_half() {};
  prob_half(float f);
  
  operator float() const;

  prob_half& operator=(prob_half h);
  prob_half& operator=(float f);

  prob_half& operator+=(prob_half h);
  prob_half& operator+=(float f);

  prob_half& operator-=(prob_half h);
  prob_half& operator-=(float f);

  prob_half& operator*=(prob_half h);
  prob_half& operator*=(float f);

  prob_half& operator/=(prob_half h);
  prob_half& operator/=(float f);

  unsigned short bits() const;
  void setBits(uint16_t bits);

  bool isNormalized() const;
  bool isDenormalized() const;
  bool isZero() const;

public:
  union float_holder {
    uint32_t i;
    float	f;
  };

  union half_holder {
    uint16_t bits;
    struct {
      unsigned m : 11;
      unsigned e : 5;
    };
  };

  static uint16_t convert(float f) {
    float_holder fh;
    fh.f = f;
    return convert(fh.i);
  }
  
private:
  static uint16_t convert(uint32_t i);

  half_holder _h;

  static const float_holder _toFloat[1 << 16];
};

#define __MASK(size) ((size) == 32 ? 0xFFFFFFFFu : (1u << (size)) - 1)
inline prob_half::prob_half(float f) {
  float_holder x;
  x.f = f;

  // Special case - zero. prob_half's do not have sign, thus we drop the sign of
  // the zero here.
  if (f == 0) {
    _h.bits = 0;
    return;
  }

  // Extract the exponent. Mantissa has 23 bits here.
  int e = ((x.i >> 23) & __MASK(8)) - 127 + 31;
  if (e > 0 && e <= 31) {
    // Simple case (normal number) - round the significand, m, to 11
    // bits and combine it with exponent.
    unsigned m = x.i & __MASK(23);
    _h.e = e;
    _h.m = (m + __MASK(23 - 11 - 1) + ((m >> (23 - 11)) & 1)) >> (23 - 11);
  } else {
    // Difficult case (denormal number) - call a function.
    _h.bits = convert(x.i);
  }
}
#undef __MASK

inline prob_half::operator float() const{
  return _toFloat[_h.bits].f;
}

inline prob_half& prob_half::operator=(prob_half h) {
  _h.bits = h._h.bits;
  return *this;
}

inline prob_half& prob_half::operator=(float f) {
  *this = prob_half(f);
  return *this;
}

inline prob_half& prob_half::operator+=(prob_half h) {
  *this = prob_half(float(*this) + float(h));
  return *this;
}

inline prob_half& prob_half::operator+=(float f) {
  *this = prob_half(float(*this) + f);
  return *this;
}

inline prob_half& prob_half::operator-=(prob_half h) {
  *this = prob_half(float(*this) - float(h));
  return *this;
}

inline prob_half& prob_half::operator-= (float f) {
  *this = prob_half(float (*this) - f);
  return *this;
}

inline prob_half &prob_half::operator*= (prob_half h) {
  *this = prob_half(float(*this) * float(h));
  return *this;
}

inline prob_half &prob_half::operator*= (float f) {
  *this = prob_half (float (*this) * f);
  return *this;
}

inline prob_half &prob_half::operator/= (prob_half h) {
  *this = prob_half (float (*this) / float (h));
  return *this;
}

inline prob_half &prob_half::operator/= (float f) {
  *this = prob_half (float (*this) / f);
  return *this;
}

inline bool prob_half::isNormalized () const {
  return _h.e > 0;
}

inline bool prob_half::isDenormalized () const {
  return _h.e == 0 && _h.m != 0;
}

inline bool prob_half::isZero () const {
  return _h.e == 0 && _h.m == 0;
}

inline uint16_t prob_half::bits() const {
  return _h.bits;
}

inline void prob_half::setBits(uint16_t bits) {
  _h.bits = bits;
}

inline std::ostream& operator<<(std::ostream &os, prob_half h) {
  os << float (h);
  return os;
}

inline std::istream& operator>>(std::istream &is, prob_half &h) {
  float f;
  is >> f;
  h = prob_half(f);
  return is;
}

#endif
