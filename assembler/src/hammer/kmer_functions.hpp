/**
 * @file    preproc.cpp
 * @author  Alex Davydow
 * @version 1.0
 *
 * @section LICENSE
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation; either version 2 of
 * the License, or (at your option) any later version.
 *
 * @section DESCRIPTION
 *
 * In this there are several functions which can be helpfull while
 * generating kmers from reads.
 *
 */
#ifndef HAMMER_KMERFUNCTIONS_HPP_
#define HAMMER_KMERFUNCTIONS_HPP_
#include <vector>
#include "common/read/read.hpp"

// Strange structure. Looks like freq is never used.
struct KMerStat {
  size_t count;
  float freq;
};
/**
 * trim bad quality nucleotides from start and end of the read
 * @return size of the read left
 */
uint32_t TrimBadQuality(Read *r, int bad_quality_threshold = 2);

/**
 * @param k k as in k-mer
 * @param start start point
 * @return the first starting point of a valid k-mer >=start; return
 * seq_.size() if no such place exists
 */
uint32_t FirstValidKmerPos(const Read &r, uint32_t start, uint32_t k);

Sequence GetSubSequence(const Read &r, uint32_t start, uint32_t length);

template<uint32_t kK>
vector< Seq<kK> > GetKMers(const Read &r) {
  const string &seq = r.getSequenceString();
  vector< Seq<kK> > ans;
  uint32_t pos = 0;
  while (true) {
    pos = FirstValidKmerPos(r, pos, kK);
    if (pos >= seq.size()) break;
    Seq<kK> kmer = Seq<kK>(GetSubSequence(r, pos, kK));
    while (true) {
      ans.push_back(kmer);
      if (pos + kK < r.size() && is_nucl(seq[pos + kK])) {
        kmer = kmer << r[pos + kK];
        ++pos;
      } else {
        pos += kK;
        break;
      }
    }
  }
  return ans;
}

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, KMerStatMap *v) {
  vector< Seq<kK> > kmers = GetKMers<kK>(r);
  for (uint32_t i = 0; i < kmers.size(); ++i) {
    ++(*v)[kmers[i]].count;
  }
}

#endif  // HAMMER_KMERFUNCTIONS_HPP_
