/***************************************************************************
 * Title:          HashUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashUtils.h"

#include "alignutils.h"
#include "SeqUtils.h"

ssize_t HashValue::size = 13;
ssize_t HashValue::indexSize = 8;
ssize_t HashValue::convertSize = 5;

ssize_t nCreated, nStored, nFilledBuckets;
Chain* HashTable::Find(HashValue hashValue, Chain* &chainPtr) {
  ssize_t bucketPos;
  chainPtr= NULL;
  assert(hashValue.GetIndex() < nBuckets);
  Chain *lastChainPtr;
  lastChainPtr = NULL;
  chainPtr = buckets[hashValue.GetIndex()];
  while (chainPtr != NULL and
	 chainPtr->hashValue != hashValue) {
    lastChainPtr = chainPtr;
    chainPtr = chainPtr->next;
  }
  if (chainPtr != NULL) {
    /*
      std::cout << "found " << chainPtr->hashValue.hashValue << " " 
      << hashValue.hashValue << " " << hashValue.GetIndex() << " " << chainPtr <<  std::endl;
    */
  }
  return lastChainPtr;
}

Chain* HashTable::StorePos(HashValue hashValue, HashPos &p) {
  assert(hashValue.GetIndex() < nBuckets);
  Chain *chainPtr, *lastChainPtr;
  chainPtr = lastChainPtr = NULL;
  ++nStored;
  lastChainPtr = Find(hashValue, chainPtr);
  if (chainPtr == NULL ) {
    // This value was not found, append it to the last chain ptr
    // should have found a last chain ptr;
    chainPtr = new Chain(hashValue);
    if (lastChainPtr == NULL) 
      buckets[hashValue.GetIndex()] = chainPtr;
    else
      lastChainPtr->next = chainPtr;
    ++nCreated;
  }
  else {
    /*
      std::cout << "appending to " << chainPtr->hashValue.hashValue << " curlen: " << chainPtr->positions.size() << " ptr: " << chainPtr << std::endl;
    */
      }
  assert(chainPtr != NULL);
  chainPtr->positions.push_back(p);
  //  std::cout << p.pos << " cp " << hashValue.hashValue << " " << chainPtr->positions.size() << std::endl;
}


ssize_t HashGenome(DNASequence &seq, char seqNumber, 
	       HashPattern &hashPattern,
	       ssize_t maskRepeats,
	       HashTable &hashTable) {
  assert(hashTable.nBuckets > 0);
  nCreated = 0; nStored = 0;
  ssize_t p;
  ssize_t lastPos = seq.length - hashPattern.mask.size();
  HashValue hashValue;
  for (p = 0; p < lastPos; p++) {
    if (GetHashValue(seq, p, hashPattern, maskRepeats, hashValue)) {
      HashPos pos(0,p);
      hashTable.StorePos(hashValue,pos);
    }
  }
  std::cout << "genome statistics: " << nCreated << " " << nStored << std::endl;
}



ssize_t GetHashValue(DNASequence &seq, ssize_t pos, 
		 HashPattern &hashPattern,
		 ssize_t maskRepeats,
		 HashValue &hashValue) {
  ssize_t i;
  ssize_t value;
  value = 0;
  // Can't mask this position
  assert(pos + hashPattern.mask.size() < seq.length);
  assert(hashPattern.maskIndices.size() > 0);
  unsigned char val;
  unsigned char nuc;
  DNASequence frag;
  char nucs[100];
  for (i = 0; i < hashPattern.maskIndices.size()-1; i++) {
    val = seq.seq[pos+hashPattern.maskIndices[i]];
    if ( nuc_index[val] > 3 or 
	 (maskRepeats and numeric_nuc_index[(unsigned char)val] < 0)) {
      return 0;
    }
    value += nuc_index[val];
    value <<= 2;
    nucs[i] = nuc_char[val];
  }
  val = seq.seq[pos+hashPattern.maskIndices[i]];
  if ( nuc_index[val] > 3 or 
       (maskRepeats and numeric_nuc_index[(unsigned char)val] < 0)) {
    return 0;
  }
  value += nuc_index[val];
  nucs[i] = nuc_char[val];
  nucs[i+1] = 0;
  hashValue = value;
  return 1;
}


void QueryHash(DNASequence &refSeq,
	       HashPattern &hashPattern,
	       HashTable &qryHashTable) {

  ssize_t pos;
  HashValue hashValue;
  ssize_t lastPos = refSeq.length - hashPattern.size();
  Chain *chain;
  for (pos = 0; pos < lastPos; pos++) {
    if (GetHashValue(refSeq, pos, hashPattern, 1, hashValue)) {
      // Try to find the reference in query
      qryHashTable.Find(hashValue, chain);
      if (chain != NULL) {
	ssize_t c;
	for (c = 0; c < chain->positions.size(); c++) {
	  std::cout << " " << chain->positions[c].pos;
	}
	std::cout << std::endl;
      }
    }
  }
}
