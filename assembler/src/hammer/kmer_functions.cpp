/*
* kmer_functions.cpp
*
*  Created on: 13.07.2011
*      Author: adavydow
*/
#include "kmer_functions.hpp"
#include <string>
#include <vector>
#include "common/read/read.hpp"
#include "common/sequence/nucl.hpp"

using std::string;
using std::vector;

Sequence GetSubSequence(const Read &r, uint32_t start, uint32_t length) {
  const string &seq = r.getSequenceString();
  assert(length > 0 && start >= 0 && start + length <= seq.size());
  if (!r.isValid()) {
    for (uint32_t i = 0; i < length; ++i) {
      assert(is_nucl(seq[start + i]));
    }
  }
  return Sequence(seq.substr(start, length));
}

uint32_t FirstValidKmerPos(const Read &r, uint32_t start, uint32_t k) {
  const string &seq = r.getSequenceString();
  uint32_t curHypothesis = start;
  uint32_t i = start;
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

uint32_t TrimBadQuality(Read *r, int bad_quality_threshold) {
  string &seq = r->seq_;
  string &qual = r->qual_;
  uint32_t start = 0;
  for (; start < qual.size(); ++start) {
    if (qual[start] > bad_quality_threshold)
      break;
  }
  seq.erase(seq.begin(), seq.begin() + start);
  qual.erase(qual.begin(), qual.begin() + start);
  if (seq.size() > 0) {
    uint32_t end = seq.size() - 1;
    for (; end > 0; --end) {
      if (qual[end] > bad_quality_threshold)
        break;
    }
    seq.erase(seq.begin() + end + 1, seq.end());
    qual.erase(qual.begin() + end + 1, qual.end());
  }
  r->valid_ = r->updateValid();
  return seq.size();
}
