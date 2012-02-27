/***************************************************************************
 * Title:          HashUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashUtils.h"

#include "SeqUtils.h"
#include "align/alignutils.h"

ssize_t HashValue::size = 13;
ssize_t HashValue::indexSize = 11;
ssize_t HashValue::convertSize = 5;

// TODO: Fix numerous compiler warnings throughout this file:
//   "warning: dereferencing pointer 'chainPtr.127' does break strict-aliasing rules"
//   "warning: dereferencing type-punned pointer will break strict-aliasing rules"

ssize_t nCreated, nStored, nFilledBuckets;
Chain* HashTable::Find(HashValue &hashValue, Chain* &chainPtr) {
  //UNUSED// ssize_t bucketPos;
  chainPtr= (Chain *) NULL;
  assert(hashValue.GetIndex() < nBuckets);
  Chain *lastChainPtr;
  lastChainPtr = (Chain *) NULL;
  chainPtr = buckets[hashValue.GetIndex()];
  while (chainPtr != (Chain *) NULL and
	 chainPtr->hashValue != hashValue) {
    lastChainPtr = chainPtr;
    chainPtr = chainPtr->next;
  }
  if (chainPtr != (Chain *) NULL) {
    /*
      std::cout << "found " << chainPtr->hashValue.hashValue << " " 
      << hashValue.hashValue << " " << hashValue.GetIndex() << " " << chainPtr <<  std::endl;
    */
  }
  return lastChainPtr;
}

//CountChain* 
void HashTable::StoreCount(HashValue &hashValue) {
  assert(hashValue.GetIndex() < nBuckets);
  CountChain *chainPtr, *lastChainPtr;
  chainPtr = lastChainPtr = (CountChain *) NULL;
  ++nStored;
  lastChainPtr = (CountChain*) Find(hashValue, (Chain*&)chainPtr);
  if (chainPtr == (CountChain *) NULL ) {
    // This value was not found, append it to the last chain ptr
    // should have found a last chain ptr;
    chainPtr = new CountChain;
    if (lastChainPtr == (CountChain *) NULL) 
      buckets[hashValue.GetIndex()] = chainPtr;
    else
      lastChainPtr->next = chainPtr;
    ++nCreated;
  }
  assert(chainPtr != (CountChain *) NULL);
  chainPtr->count++;
}


//PositionChain*
void HashTable::StorePos(HashValue &hashValue, HashPos &p) {
  assert(hashValue.GetIndex() < nBuckets);
  PositionChain *chainPtr, *lastChainPtr;
  chainPtr = lastChainPtr = (PositionChain *) NULL;
  ++nStored;
  //UNUSED// Chain *tmpPtr;
  lastChainPtr = (PositionChain*) Find(hashValue, (Chain*&) chainPtr);
  if (chainPtr == (PositionChain *) NULL ) {
    // This value was not found, append it to the last chain ptr
    // should have found a last chain ptr;
    chainPtr = new PositionChain(hashValue);
    if (lastChainPtr == (PositionChain *) NULL) 
      buckets[hashValue.GetIndex()] = chainPtr;
    else
      lastChainPtr->next = chainPtr;
    ++nCreated;
  }
  assert(chainPtr != (PositionChain *) NULL);
  chainPtr->positions.push_back(p);
}

ssize_t CountHashGenome(DNASequence &seq, char seqNumber, 
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
      hashTable.StoreCount(hashValue);
    }
  }
  return 1;
}

void HashGenome(DNASequence &seq, char seqNumber, 
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
  size_t value;
  value = 0;
  // Can't mask this position
  assert(pos + hashPattern.mask.size() <= seq.length);
  assert(hashPattern.maskIndices.size() > 0);
  unsigned char val;
  unsigned char nuc;
	char maskNuc;
  DNASequence frag;

	// TODO: has overlap with IntegralTuple.h
  for (i = 0; i < hashPattern.maskIndices.size(); i++) {
    val = seq.seq[pos+hashPattern.maskIndices[i]];
		nuc = numeric_nuc_index[(unsigned char) val];
		maskNuc = unmasked_nuc_index[(unsigned char) val];
    if ( (!maskRepeats and nuc > 3) or 
				 (maskRepeats and maskNuc < 0)) {
      return 0;
    }
    value <<= 2;
    value += nuc;
  }
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
    }
  }
}
