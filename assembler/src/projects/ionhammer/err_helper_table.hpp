//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_ERR_HELPER_TABLE_HPP__
#define __HAMMER_ERR_HELPER_TABLE_HPP__

#include "hkmer.hpp"

#include <cassert>
#include <cstdlib>
#include <istream>
#include <string>
#include <vector>

#include "utils/logger/logger.hpp"

namespace hammer {

namespace errHelper {

/// Type of error
enum Hint { kMismatch, kInsertion, kDeletion };

namespace internal {

// maximum size of K-mers in the helper tables
static const unsigned int MAX_K = 5;

struct HelperTable {
  const unsigned k_;
  const uint32_t *storage_;

  template <typename It1, typename It2>
  Hint lookupHint(const It1 &x_it, const It2 &y_it, size_t x_nfront,
                  size_t y_nfront) const {
    VERIFY(k_ <= MAX_K);
    unsigned x_code = getCode(x_it, x_nfront, k_);
    unsigned y_code = getCode(y_it, y_nfront, k_);

    unsigned code = x_code + (y_code << (2 * k_));
    uint32_t bt = storage_[code / 16];  // 16 hints per uint32_t
    unsigned shift = (code % 16) * 2;
    return static_cast<Hint>((bt >> shift) & 0x3);
  }

  template <typename HRunIter>
  static unsigned getCode(const HRunIter &x_it, size_t x_nfront, size_t k) {
    unsigned code = 0;
    unsigned len = 0;
    auto nucl = x_it->nucl;
    for (len = 0; len < x_nfront && len < k; ++len) code |= nucl << (2 * len);

    if (len == k) return code;

    for (HRunIter it = x_it + 1;; ++it) {
      for (size_t i = 0; i < it->len; ++i) {
        code |= it->nucl << (2 * len++);
        if (len == k) return code;
      }
    }

    assert(false);
  }
};

// tables for k = 1, 2, ..., MAX_K
extern const HelperTable helper_tables[];

template <typename HRunIter>
static inline size_t getNumberOfRemainingBases(const HRunIter &x_it,
                                               const HRunIter &x_end,
                                               size_t x_nfront) {
  size_t n = x_nfront;
  if (n >= MAX_K) return MAX_K;

  for (HRunIter it = x_it + 1; it != x_end; ++it) {
    n += it->len;
    if (n >= MAX_K) return MAX_K;
  }

  return n;
}

};  // namespace internal

/// Estimate what kind of error occurred at the position
template <typename It1, typename It2>
static inline Hint getHint(const It1 &x_begin, const It1 &x_end,
                           const It2 &y_begin, const It2 &y_end,
                           size_t x_nfront, size_t y_nfront) {
  VERIFY(x_nfront <= x_begin->len);
  VERIFY(y_nfront <= y_begin->len);
  size_t x_rem = internal::getNumberOfRemainingBases(x_begin, x_end, x_nfront);
  size_t y_rem = internal::getNumberOfRemainingBases(y_begin, y_end, y_nfront);

  auto &table = internal::helper_tables[std::min(x_rem, y_rem) - 1];
  return table.lookupHint<It1, It2>(x_begin, y_begin, x_nfront, y_nfront);
}

};  // namespace errHelper
};  // namespace hammer

#endif  // __HAMMER_ERR_HELPER_TABLE_HPP__
