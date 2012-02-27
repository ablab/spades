/***************************************************************************
 * Title:          HashUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef HASHUTILS_H_
#define HASHUTILS_H_

#include "mctypes.h"
#include "DNASequence.h"
#include <vector>


class HashValue {
public:
  static ssize_t size;
  static ssize_t indexSize;
  static ssize_t convertSize;
  size_t hashValue;
  ssize_t GetIndex() {
    return hashValue >> (convertSize*2);
  }
  HashValue & operator=(const HashValue &hv) {
		if (this != &hv) {
			hashValue = hv.hashValue;
		}
    return *this;
  }
  HashValue & operator=(const ssize_t &val) {
		hashValue = val;
    return *this;
  }
  bool operator==(const HashValue &val) const {
    return val.hashValue == hashValue;
  }
  bool operator!=(const HashValue &val) const {
    return val.hashValue != hashValue;
  }
};

class HashPattern {
public:
  BitVector mask;
  IntVector maskIndices;
  ssize_t nChars;
  ssize_t length;
  void SetNChars() {
    ssize_t i;
    nChars = 0;
    for (i = 0; i < mask.size(); i++) {
      if (mask[i]) ++nChars;
    }
    HashValue::size = 2*nChars;
    HashValue::convertSize = HashValue::size - HashValue::indexSize;
    if (HashValue::convertSize < 0) 
      HashValue::convertSize = 0;
  }
  void SetUngappedSize(ssize_t size) {
    ssize_t i;
    mask.clear();
    maskIndices.clear();
    for (i = 0; i < size; i++ ) {
      mask.push_back(1);
      maskIndices.push_back(i);
    }
    length = maskIndices.size();
  }
  ssize_t size() { return mask.size(); }
};

class HashPos {
public:
  char sequence;
  ssize_t pos;
  HashPos() {
    sequence = -1;
    pos = -1;
  }
  HashPos(char s, ssize_t p) {
    sequence = s;
    pos = p;
  }
  HashPos &operator=(const HashPos &p2) {
		if (this != &p2) {
			sequence = p2.sequence;
			pos = p2.pos;
		}
    return *this;
  }
};


class Chain {
public:
  HashValue hashValue;
  Chain *next;
  Chain(HashValue &hv) {
    hashValue = hv;
    next = NULL;
  }
  Chain() {
    next = NULL;
    hashValue = -1;
  }
};

class PositionChain : public Chain {
public:
  std::vector<HashPos> positions;
  PositionChain(HashValue &hv) : Chain(hv) {};
  PositionChain() : Chain() {};
};


class CountChain : public Chain {
public:
  ssize_t count;
  CountChain() : Chain(), count(0) {}
};
  
class HashTable {
public:
  ssize_t nBuckets;
  std::vector<Chain*> buckets;
  HashTable(ssize_t nBuckets) {
    nBuckets = CalcNBuckets(nBuckets);
    Init(nBuckets);
  }
  HashTable() {
    nBuckets = CalcNBuckets(HashValue::indexSize);
    Init(nBuckets);
  }
  void Init(ssize_t nb) {
    nBuckets = nb;
    buckets.resize(nBuckets);
    ssize_t i;
    for (i = 0; i < nBuckets; i++) {
      buckets[i] = NULL;
    }
  }
  static ssize_t CalcNBuckets(ssize_t nChars) {
    return 1 << (2*nChars);
  }

	//  PositionChain* StorePos(HashValue &hashValue, HashPos &p);
	void StorePos(HashValue &hashValue, HashPos &p);
	//  CountChain* StoreCount(HashValue &hashValue);
	void StoreCount(HashValue &hashValue);
  Chain* Find(HashValue &hashValue, Chain* &chainPtr);
};


ssize_t GetHashValue(DNASequence &seq, ssize_t pos, 
			  HashPattern &hashPattern,
			  ssize_t maskRepeats,
			  HashValue &hashValue);

ssize_t CountHashGenome(DNASequence &seq, char seqNumber, 
		    HashPattern &hashPattern,
		    ssize_t maskRepeats,
		    HashTable &hashTable);

void HashGenome(DNASequence &seq, char seqNumber, 
	       HashPattern &hashPattern,
	       ssize_t maskRepeats,
	       HashTable &hashTable);

void QueryHash(DNASequence &refSeq,
	       HashPattern &hashPattern,
	       HashTable &qryHashTable);
	       


#endif
