//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef HAMMER_KMERFREQINFO_HPP_
#define HAMMER_KMERFREQINFO_HPP_
#include <string>
struct KMerFreqInfo {
  std::string kmer;
  uint32_t count;
  float q_count;
  KMerFreqInfo() : kmer(""), count(0), q_count(0) {}
};

#endif //  HAMMER_KMERFREQINFO_HPP_
