//===- llvm/ADT/PointerEmbeddedInt.h ----------------------------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//

#ifndef LLVM_ADT_POINTEREMBEDDEDINT_H
#define LLVM_ADT_POINTEREMBEDDEDINT_H

#include "PointerLikeTypeTraits.h"
#include <climits>

namespace llvm {

/// Utility to embed an integer into a pointer-like type. This is specifically
/// intended to allow embedding integers where fewer bits are required than
/// exist in a pointer, and the integer can participate in abstractions along
/// side other pointer-like types. For example it can be placed into a \c
/// PointerSumType or \c PointerUnion.
///
/// Note that much like pointers, an integer value of zero has special utility
/// due to boolean conversions. For example, a non-null value can be tested for
/// in the above abstractions without testing the particular active member.
/// Also, the default constructed value zero initializes the integer.
template <typename IntT, int Bits = sizeof(IntT) * CHAR_BIT>
class PointerEmbeddedInt {
  uintptr_t Value;

  static_assert(Bits < sizeof(uintptr_t) * CHAR_BIT,
                "Cannot embed more bits than we have in a pointer!");

  enum : uintptr_t {
    // We shift as many zeros into the value as we can while preserving the
    // number of bits desired for the integer.
    Shift = sizeof(uintptr_t) * CHAR_BIT - Bits,

    // We also want to be able to mask out the preserved bits for asserts.
    Mask = static_cast<uintptr_t>(-1) << Bits
  };

  friend class PointerLikeTypeTraits<PointerEmbeddedInt>;

  explicit PointerEmbeddedInt(uintptr_t Value) : Value(Value) {}

public:
  PointerEmbeddedInt() : Value(0) {}

  PointerEmbeddedInt(IntT I) : Value(static_cast<uintptr_t>(I) << Shift) {
    assert((I & Mask) == 0 && "Integer has bits outside those preserved!");
  }

  PointerEmbeddedInt &operator=(IntT I) {
    assert((I & Mask) == 0 && "Integer has bits outside those preserved!");
    Value = static_cast<uintptr_t>(I) << Shift;
  }

  // Note that this imilict conversion additionally allows all of the basic
  // comparison operators to work transparently, etc.
  operator IntT() const { return static_cast<IntT>(Value >> Shift); }
};

// Provide pointer like traits to support use with pointer unions and sum
// types.
template <typename IntT, int Bits>
class PointerLikeTypeTraits<PointerEmbeddedInt<IntT, Bits>> {
  typedef PointerEmbeddedInt<IntT, Bits> T;

public:
  static inline void *getAsVoidPointer(const T &P) {
    return reinterpret_cast<void *>(P.Value);
  }
  static inline T getFromVoidPointer(void *P) {
    return T(reinterpret_cast<uintptr_t>(P));
  }
  static inline T getFromVoidPointer(const void *P) {
    return T(reinterpret_cast<uintptr_t>(P));
  }

  enum { NumLowBitsAvailable = T::Shift };
};

}

#endif
