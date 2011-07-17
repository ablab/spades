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

int32_t getKmerAnew( const std::string & seq, int32_t pos, KMer & kmer ) {
	int32_t curHypothesis = pos;
	int32_t i = pos;

	//cout << "getAnew: " << seq.data() << " " << i << endl;
	for (; i < seq.size(); ++i) {
		if (i >= K + curHypothesis) {
			kmer = KMer( seq.data() + curHypothesis, false );
			return curHypothesis + 1;
		}
		if (!is_nucl(seq[i])) {
			curHypothesis = i + 1;
		}
	}
	if (i >= K + curHypothesis) {
		kmer = KMer( seq.data() + curHypothesis, false );
		return curHypothesis + 1;
	}
	return -1;
}

int32_t NextValidKmer(const Read &r, int32_t pos, KMer & kmer) {	
	const std::string &seq = r.getSequenceString();
	if (pos == -1) { // need to get first valid kmer
		return getKmerAnew(seq, 0, kmer);
	} else {
		if (pos + K < r.size() && is_nucl(seq[pos + K])) {
			kmer = kmer << r[pos + K];
			return (pos+1);
		} else {
			return getKmerAnew(seq, pos, kmer);
		}
	}
}

vector<KMer> GetKMers(const Read &r) {
  const string &seq = r.getSequenceString();
  vector<KMer> ans;
  uint32_t pos = 0;
  while (true) {
    pos = FirstValidKmerPos(r, pos, K);
    if (pos >= seq.size()) break;
    KMer kmer = KMer(GetSubSequence(r, pos, K));
    while (true) {
      ans.push_back(kmer);
      if (pos + K < r.size() && is_nucl(seq[pos + K])) {
        kmer = kmer << r[pos + K];
        ++pos;
      } else {
        pos += K;
        break;
      }
    }
  }
  return ans;
}

/**
 * add k-mers from read to map
 */
void AddKMers(const Read &r, uint64_t readno, KMerStatMap *v) {
	int32_t pos = -1;
	KMer kmer;
	while ( (pos = NextValidKmer(r, pos, kmer)) >= 0 ) {
		++(*v)[kmer].count;
		(*v)[kmer].pos.push_back( make_pair(readno, pos - 1) );
	}
}


