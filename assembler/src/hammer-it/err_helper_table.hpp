#ifndef __HAMMER_ERR_HELPER_TABLE_HPP__
#define __HAMMER_ERR_HELPER_TABLE_HPP__

#include "hkmer.hpp"

#include <vector>
#include <istream>
#include <string>
#include <cstdlib>
#include <cassert>

#include "logger/logger.hpp"

namespace hammer {

namespace errHelper {

/// Type of error
enum Hint {
  kMismatch,
  kInsertion,
  kDeletion
};

namespace internal {

// maximum size of K-mers in the helper tables
static const unsigned int MAX_K = 5;

struct HelperTable {
  const unsigned k_;
  const uint32_t* storage_;

  template <typename T1, typename T2>
  Hint lookupHint(const T1 &x, const T1 &y,
                  size_t x_offset, size_t y_offset,
                  size_t x_nfront, size_t y_nfront) const {

    VERIFY(k_ <= MAX_K);
    unsigned x_code = getCode(x, x_offset, x_nfront, k_);
    unsigned y_code = getCode(y, y_offset, y_nfront, k_);

    unsigned code = x_code + (y_code << (2 * k_));
    uint32_t bt = storage_[code / 16]; // 16 hints per uint32_t
    unsigned shift = (code % 16) * 2;
    return static_cast<Hint>((bt >> shift) & 0x3);
  }

  template <typename Runs>
  static unsigned getCode(const Runs& x, size_t x_offset,
                          size_t x_nfront, size_t k) {
    unsigned code = 0;
    unsigned len = 0;
    auto nucl = x[x_offset].nucl;
    for (len = 0; len < x_nfront && len < k; ++len)
      code |= nucl << (2 * len);

    if (len == k)
      return code;

    for (size_t offset = x_offset + 1; ; ++offset) {
      for (size_t i = 0; i < x[offset].len; ++i) {
        code |= x[offset].nucl << (2 * len++);
        if (len == k)
          return code;
      }
    }

    assert(false);
  }
};

// tables for k = 1, 2, ..., MAX_K
extern const HelperTable helper_tables[];

template <typename Runs>
static inline size_t getNumberOfRemainingBases(const Runs &x,
                                               size_t x_offset,
                                               size_t x_nfront) {
  size_t n = x_nfront;
  if (n >= MAX_K)
    return MAX_K;

  auto sz = hammer::internal::getSize(x);
  for (size_t i = x_offset + 1; i < sz; ++i) {
    n += x[i].len;
    if (n >= MAX_K)
      return MAX_K;
  }

  return n;
}

}; // namespace internal

/// Estimate what kind of error occurred at the position
template <typename T1, typename T2>
static inline Hint getHint(const T1 &x, const T2 &y,
                           size_t x_offset, size_t y_offset,
                           size_t x_nfront, size_t y_nfront) {
  size_t x_rem = internal::getNumberOfRemainingBases(x, x_offset, x_nfront);
  size_t y_rem = internal::getNumberOfRemainingBases(y, y_offset, y_nfront);

  auto& table = internal::helper_tables[std::min(x_rem, y_rem) - 1];
  return table.lookupHint<T1, T2>(x, y, x_offset, y_offset, x_nfront, y_nfront);
}

}; // namespace errHelper
}; // namespace hammer

#endif // __HAMMER_ERR_HELPER_TABLE_HPP__
