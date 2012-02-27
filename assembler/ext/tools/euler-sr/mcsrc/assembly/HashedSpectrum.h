/***************************************************************************
 * Title:          HashedSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef HASHED_SPECTRUM_H_
#define HASHED_SPECTRUM_H_

#include "Spectrum.h"
#include "ReadPos.h"
#include "SimpleSequence.h"
#include "PagedList.h"
#include "hash/PagedHashTable.h"

class HashValueFunctor {
public:
  ssize_t hashLength;
  SimpleSequenceList *sequences;

  ssize_t MaxHashValue() {
    return 1 << (hashLength*2);
  }
	HashValueFunctor &operator=(const HashValueFunctor &hashFunction) {
		if (this != &hashFunction) {
			hashLength = hashFunction.hashLength;
			sequences  = hashFunction.sequences;
		}
		return *this;
	}
  ssize_t operator()(CountedReadPos &p) {
    ssize_t hashValue = 0;
    ssize_t i;
    ssize_t nucValue;
    for (i = 0; i < hashLength; i++ ) {
      nucValue = unmasked_nuc_index[(*sequences)[p.read].seq[p.pos + i]];
      if ( nucValue < 4) {
				hashValue <<=2;
				hashValue += nucValue;
      }
      else {
				return -1;
      }
    }
    return hashValue;
  }
};


typedef PagedHashTable<CountedReadPos, 15, HashValueFunctor, UpdateFunctor> ReadPosHashTable;


class HashedSpectrum  {
 public:

	ssize_t hashLength;
	HashedSpectrum(HashValueFunctor &calcHashValue);
	HashValueFunctor hashFunction;
	ReadPosHashTable hashTable;
	void SetHashSize(ssize_t size) {
		HashedSpectrum::hashLength = size;
	}
	
	
	ssize_t HashSequence(SimpleSequenceList &seqListPtr, ssize_t curSeq,
									 ssize_t trimFront = 0,
									 ssize_t trimEnd   = 0);

	void StoreSpectrum(SimpleSequenceList &seqList, 
										 ssize_t trimFront = 0, ssize_t trimEnd   = 0);

	void HashToReadPositionList(SimpleSequenceList &sequences,
															ReadPositions &readPos );
};


#endif
