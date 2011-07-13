/*
* kmer_functions.cpp
*
*  Created on: 13.07.2011
*      Author: adavydow
*/

#include "kmer_functions.hpp"
#include <string>
#include <utility>
#include "common/sequence/nucl.hpp"

Sequence GetSubSequence(const Read &r, size_t start, size_t length) {
  const std::string &seq = r.getSequenceString();
  assert(length > 0 && start >= 0 && start + length <= seq.size());
  if (!r.isValid()) {
    for (size_t i = 0; i < length; ++i) {
      assert(is_nucl(seq[start + i]));
    }
  }
  return Sequence(seq.substr(start, length));
}

size_t FirstValidKmerPos(const Read &r, size_t start, size_t k) {
  const std::string &seq = r.getSequenceString();
  size_t curHypothesis = start;
  size_t i = start;
  for (; i < seq.size(); ++i) {
    if (i >= k + curHypothesis)
      return curHypothesis;
    if (!is_nucl(seq[i])) {
      curHypothesis = i + 1;
    }
  }
  if (i >= k + curHypothesis) {
    return curHypothesis;
  }
  return seq.size();
}

size_t TrimBadQuality(Read &r, int bad_quality_threshold) {
  std::string &seq = r.seq_;
  std::string &qual = r.qual_;
  size_t start = 0;
  for (; start < qual.size(); ++start) {
    if (qual[start] > bad_quality_threshold)
      break;
  }
  seq.erase(seq.begin(), seq.begin() + start);
  qual.erase(qual.begin(), qual.begin() + start);
  size_t end = seq.size() - 1;
  for (; end > 0; --end) {
    if (qual[end] > bad_quality_threshold)
      break;
  }
  seq.erase(seq.begin() + end + 1, seq.end());
  qual.erase(qual.begin() + end + 1, qual.end());
  r.valid_ = r.updateValid();
  return seq.size();
}

/**
 * add k-mers from read to map
 */
void AddKMers(const Read &r, KMerStatMap &v) {
  KMerStatMap::iterator it;
  const std::string &seq = r.getSequenceString();
  size_t pos = 0;
  while (true) {
    pos = FirstValidKmerPos(r, pos, K);
    if (pos >= seq.size()) break;
    KMer kmer = KMer(GetSubSequence(r, pos, K));
    while (true) {
      ++v[kmer].count;
      if (pos + K < (int)r.size() && is_nucl(seq[pos + K])) {
	kmer = kmer << r[pos + K];
	++pos;
      } else {
	pos += K;
	break;
      }
    }
  }
}

