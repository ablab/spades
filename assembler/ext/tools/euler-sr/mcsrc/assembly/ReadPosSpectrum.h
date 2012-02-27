/***************************************************************************
 * Title:          ReadPosSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/05/2009
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

class HashValueFunctor {
public:
  ssize_t hashLength;
  SimpleSequenceList *sequences;

  ssize_t MaxHashValue() {
    return 1 << (hashLength*2);
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


class ReadPosSpectrum : public Spectrum {
 public:
	static ssize_t hashLength;
	HashValueFunctor hash;
	ReadPosHashTable hashTable;
	void SetHashSize(ssize_t size) {
		HashSpectrum::hashLength = size;
		hash.hashLength = size;
	}
	
	void Read(std::string &fileName, ssize_t minMult = 0);
	ssize_t Write(std::string &fileName, ssize_t minMult = 0);
	ssize_t IncrementKmer(T &tuple);
	ssize_t FindKmer(T &tuple);

}

ssize_t HashSequence(SimpleSequenceList &seqListPtr, ssize_t curSeq,
								 ReadPosHashTable &hashTable, HashValueFunctor calcHashValue,
								 ssize_t trimFront = 0,
								 ssize_t trimEnd   = 0);

void StoreSpectrum(SimpleSequenceList &seqList, 
									 ReadPosHashTable &hashTable, 
									 HashValueFunctor &calcHashValue,
									 ssize_t trimFront = 0,
									 ssize_t trimEnd   = 0);

void HashToReadPositionList(ReadPosHashTable &hashTable,
														SimpleSequenceList &sequences,
														ssize_t hashLength,
														ReadPositions &readPos );


#endif
