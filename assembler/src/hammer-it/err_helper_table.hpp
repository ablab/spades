#ifndef __HAMMER_ERR_HELPER_TABLE_HPP__
#define __HAMMER_ERR_HELPER_TABLE_HPP__

#include "hkmer.hpp"

#include <vector>
#include <istream>
#include <string>
#include <cstdlib>
#include <cassert>

namespace hammer {

namespace errHelper {

/// Type of error
enum Hint {
  kMismatch,
  kInsertion,
  kDeletion
};

namespace internal {
class HelperTable {
 public:
  // Load table from file stream
  HelperTable(unsigned k, std::istream& stream);

  Hint lookupHint(const hammer::HKMer& x, const hammer::HKMer& y,
                  size_t x_offset, size_t y_offset,
                  size_t x_nfront, size_t y_nfront) const {

    unsigned x_code = getCode(x, x_offset, x_nfront, k_);
    unsigned y_code = getCode(y, y_offset, y_nfront, k_);

    unsigned code = x_code + (y_code << (2 * k_));
    char bt = storage_[code / 4];
    unsigned shift = (code % 4) * 2;
    return static_cast<Hint>((bt >> shift) & 0x3);
  }

 private:
  unsigned k_;
  std::vector<char> storage_;

  static unsigned getCode(const hammer::HKMer& x, size_t x_offset,
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

// maximum size of K-mers in the helper tables
static const unsigned int MAX_K = 5;

// tables for k = 1, 2, ..., MAX_K
extern std::vector<HelperTable> helper_tables;

static inline size_t getNumberOfRemainingBases(const hammer::HKMer& x,
                                               size_t x_offset,
                                               size_t x_nfront) {
  size_t n = x_nfront;
  if (n >= MAX_K)
    return MAX_K;

  for (size_t i = x_offset + 1; i < hammer::K; ++i) {
    n += x[i].len;
    if (n >= MAX_K)
      return MAX_K;
  }

  return n;
}

}; // namespace internal

/// Load tables from file
void initHelperTables(const std::string& filename);

/// Estimate what kind of error occurred at the position
static inline Hint getHint(const hammer::HKMer& x, const hammer::HKMer& y,
                           size_t x_offset, size_t y_offset,
                           size_t x_nfront, size_t y_nfront) {
  size_t x_rem = internal::getNumberOfRemainingBases(x, x_offset, x_nfront);
  size_t y_rem = internal::getNumberOfRemainingBases(y, y_offset, y_nfront);

  auto& table = internal::helper_tables[std::min(x_rem, y_rem) - 1];
  return table.lookupHint(x, y, x_offset, y_offset, x_nfront, y_nfront);
}

}; // namespace errHelper
}; // namespace hammer

#endif // __HAMMER_ERR_HELPER_TABLE_HPP__
