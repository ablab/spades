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
int32_t getKmerAnew( const std::string & seq, int32_t pos, Seq<kK> & kmer ) {
  int32_t curHypothesis = pos;
  int32_t i = pos;
  for (; i < seq.size(); ++i) {
    if (i >= kK + curHypothesis) {
      kmer = Seq<kK>(seq.data() + curHypothesis, false);
      return curHypothesis + 1;
    }
    if (!is_nucl(seq[i])) {
      curHypothesis = i + 1;
    }
  }
  if (i >= kK + curHypothesis) {
    kmer = Seq<kK>(seq.data() + curHypothesis, false);
    return curHypothesis + 1;
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
int32_t NextValidKmer(const Read &r, int32_t pos, Seq<kK> & kmer) {
  const std::string &seq = r.getSequenceString();
  if (pos == -1) { // need to get first valid kmer
    return getKmerAnew<kK>(seq, 0, kmer);
  } else {
    if (pos + kK < r.size() && is_nucl(seq[pos + kK])) {
      kmer = kmer << r[pos + kK];
      return (pos + 1);
    } else {
      return getKmerAnew<kK>(seq, pos, kmer);
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
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v) {
  int32_t pos = -1;
  Seq<kK> kmer;
  while ( (pos = NextValidKmer<kK>(r, pos, kmer)) >= 0 ) {
    ++(*v)[kmer].count;
    (*v)[kmer].pos.push_back( make_pair(readno, pos - 1) );
  }  
}
#endif  // HAMMER_KMERFUNCTIONS_HPP_
