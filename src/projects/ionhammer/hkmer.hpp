//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_HKMER_HPP__
#define __HAMMER_HKMER_HPP__

#include <cstdlib>
#include <limits>
#include "HSeq.hpp"

namespace hammer {

const uint32_t K = 16;
using HKMer = HSeq<K>;

struct HKMerDistanceResult {
  double hamming_ = 0;
  double levenshtein_ = 0;

  HKMerDistanceResult(double hamming = 0, double lev = 0)
      : hamming_(hamming), levenshtein_(lev) {}
};

inline HKMerDistanceResult hkmerDistance(const HKMer& left,
                                         const HKMer& right) {
  HKMerDistanceResult dist = {0, 0};

  for (uint32_t i = 0; i < K; ++i) {
    if (left[i].nucl != right[i].nucl) {
      return {std::numeric_limits<double>::infinity(),
              std::numeric_limits<double>::infinity()};
    }

    if (left[i].len != right[i].len) {
      dist.hamming_ += 1;
      dist.levenshtein_ += std::abs(left[i].len - right[i].len);
    }
  }
  return dist;
}



};  // namespace hammer

#endif  // __HAMMER_HKMER_HPP__
