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
#include "hammer/kmer_stat.hpp"

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

/**
 * get next valid kmer in a new position
 */
template<uint32_t kK>
int32_t getKmerAnew(const Read &r, int32_t pos, Seq<kK> & kmer ) {
  const std::string &seq = r.getSequenceString();
  int32_t curHypothesis = pos;
  int32_t i = pos;
  for (; i < seq.size(); ++i) {
    if (i >= kK + curHypothesis) {
      break;
    }
    if (!is_nucl(seq[i])) {
      curHypothesis = i + 1;
    }
  }
  if (i >= kK + curHypothesis) {
    kmer = Seq<kK>(seq.data() + curHypothesis, false);
    return curHypothesis;
  }
  return -1;
}

/**
 * @param kmer get next valid k-mer
 * @param pos starting point
 * @return the first starting point of a valid k-mer >=start; return
 * -1 if no such place exists
 */
template<uint32_t kK>
int32_t NextValidKmer(const Read &r, int32_t prev_pos, Seq<kK> & kmer) {
  const std::string &seq = r.getSequenceString();
  // cout << "  run NextValidKmer(" << seq.data() << ", " << prev_pos << ")" << endl;
  if (prev_pos == -1) { // need to get first valid kmer
    return getKmerAnew<kK>(r, 0, kmer);
  } else {
    if (prev_pos + kK < r.size() && is_nucl(seq[prev_pos + kK])) {
      kmer = kmer << r[prev_pos + kK];
      return (prev_pos + 1);
    } else {
      return getKmerAnew<kK>(r, prev_pos + 1, kmer);
    }
  }
}

template<uint32_t kK>
vector< Seq<kK> > GetKMers(const Read &r) {
  vector< Seq<kK> > ret;
  int32_t pos = -1;
  Seq<kK> kmer;
  while ((pos = NextValidKmer<kK>(r, pos, kmer)) >= 0) {
    ret.push_back(kmer);
  }  
  return ret;
}

/**
 * add k-mers from read to map
 */
template<uint32_t kK, typename KMerStatMap>
void AddKMers(const Read &r, int64_t readno, KMerStatMap *v) {
  int32_t pos = -1;
  Seq<kK> kmer;
  while ( (pos = NextValidKmer<kK>(r, pos, kmer)) >= 0 ) {
    // cout << pos << "  " << kmer.str().data() << endl;
    ++(*v)[kmer].count;
    (*v)[kmer].pos.push_back( make_pair(readno, pos) );
  }  
}
#endif  // HAMMER_KMERFUNCTIONS_HPP_
