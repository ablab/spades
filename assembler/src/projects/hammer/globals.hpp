//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_GLOBALS_HPP_
#define HAMMER_GLOBALS_HPP_

#include "kmer_stat.hpp"

class KMerData;

struct Globals {
  static int iteration_no;

  static std::vector<uint32_t> * subKMerPositions;
  static KMerData *kmer_data;

  static char char_offset;
  static bool char_offset_user;

  static double quality_probs[256];
  static double quality_lprobs[256];
  static double quality_rprobs[256];
  static double quality_lrprobs[256];
};

inline double getProb(const KMerStat &kmc, size_t i, bool log) {
  uint8_t qual = getQual(kmc, i);

  return (log ? Globals::quality_lprobs[qual] : Globals::quality_probs[qual]);
}

inline double getRevProb(const KMerStat &kmc, size_t i, bool log) {
  uint8_t qual = getQual(kmc, i);

  return (log ? Globals::quality_lrprobs[qual] : Globals::quality_rprobs[qual]);
}

#endif //  HAMMER_GLOBALS_HPP_

