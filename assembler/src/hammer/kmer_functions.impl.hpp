/**
 * @file    kmer_functions.impl.hpp
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
#include <string>

/**
 * get next valid kmer in a new position
 */
template<uint32_t kK>
int getKmerAnew(const Read &r, uint32_t pos, Seq<kK> *kmer);

/**
 * get next valid kmer in a new position
 */
template<uint32_t kK>
int getKmerAnew(const Read &r, uint32_t pos, Seq<kK> *kmer) {
  const std::string &seq = r.getSequenceString();
  uint32_t curHypothesis = pos;
  uint32_t i = pos;
  for (; i < seq.size(); ++i) {
    if (i >= kK + curHypothesis) {
      break;
    }
    if (!is_nucl(seq[i])) {
      curHypothesis = i + 1;
    }
  }
  if (i >= kK + curHypothesis) {
    *kmer = Seq<kK>(seq.data() + curHypothesis, false);
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
int NextValidKmer(const Read &r, int prev_pos, Seq<kK> *kmer) {
  const std::string &seq = r.getSequenceString();
  if (prev_pos == -1) { // need to get first valid kmer
    return getKmerAnew<kK>(r, 0, kmer);
  } else {
    if (prev_pos + kK < r.size() && is_nucl(seq[prev_pos + kK])) {
      *kmer = *kmer << r[prev_pos + kK];
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
  while ((pos = NextValidKmer<kK>(r, pos, &kmer)) >= 0) {
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
  while ( (pos = NextValidKmer<kK>(r, pos, &kmer)) >= 0 ) {
    ++(*v)[kmer].count;
    (*v)[kmer].pos.push_back( make_pair(readno, pos - 1) );
  }  
}
