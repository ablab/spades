//===-- Atomic.cpp - Atomic Operations --------------------------*- C++ -*-===//
//
//                     The LLVM Compiler Infrastructure
//
// This file is distributed under the University of Illinois Open Source
// License. See LICENSE.TXT for details.
//
//===----------------------------------------------------------------------===//
//
//  This file implements atomic operations.
//
//===----------------------------------------------------------------------===//

#include "llvm/Support/Atomic.h"

using namespace llvm;

void sys::MemoryFence() {
  __sync_synchronize();
}

sys::cas_flag sys::CompareAndSwap(volatile sys::cas_flag* ptr,
                                  sys::cas_flag new_value,
                                  sys::cas_flag old_value) {
  return __sync_val_compare_and_swap(ptr, old_value, new_value);
}

sys::cas_flag sys::AtomicIncrement(volatile sys::cas_flag* ptr) {
  return __sync_add_and_fetch(ptr, 1);
}

sys::cas_flag sys::AtomicDecrement(volatile sys::cas_flag* ptr) {
  return __sync_sub_and_fetch(ptr, 1);
}

sys::cas_flag sys::AtomicAdd(volatile sys::cas_flag* ptr, sys::cas_flag val) {
  return __sync_add_and_fetch(ptr, val);
}

sys::cas_flag sys::AtomicMul(volatile sys::cas_flag* ptr, sys::cas_flag val) {
  sys::cas_flag original, result;
  do {
    original = *ptr;
    result = original * val;
  } while (sys::CompareAndSwap(ptr, result, original) != original);

  return result;
}

sys::cas_flag sys::AtomicDiv(volatile sys::cas_flag* ptr, sys::cas_flag val) {
  sys::cas_flag original, result;
  do {
    original = *ptr;
    result = original / val;
  } while (sys::CompareAndSwap(ptr, result, original) != original);

  return result;
}
