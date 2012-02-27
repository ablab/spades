/***************************************************************************
 * Title:          ReadPosSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashedSpectrum.h"


ssize_t CountedReadPos::hashLength;
SimpleSequenceList* CountedReadPos::sequences;

template<typename T>
void Read(std::string &fileName, ssize_t minMult) {
	std::cout << "Reading a hash is not supported" << std::endl;
	exit(1);
	// return 0;
}


template<typename T>
ssize_t Write(std::string &fileName, ssize_t minMult) {
	std::cout << "Writing a hash is not supported." << std::endl;
	exit(1);
	return 0;
}


template<typename T>
ssize_t HashedSpectrum<T>::IncrementKmer(T &tuple) {
	if (tuple.Valid()) {
		if (hashTable.Store(tuple, calHashValue));
	}
}


template<typename T>
ssize_t HashedSpectrum<T>::FindKmer(T &tuple) {
	if (hashTable.find(Tuple))
		return 1;
}


ssize_t HashSequence(SimpleSequenceList &seqListPtr, ssize_t curSeq,
								 ReadPosHashTable &hashTable, HashValueFunctor calcHashValue, 
								 ssize_t trimFront, ssize_t trimEnd) {
  ssize_t p, p2;
  CountedReadPos readPos;
  ssize_t newTupleStored = 0;
  ssize_t valid = 1;
  readPos.read = curSeq;
  for (p = trimFront; p < seqListPtr[curSeq].length - CountedReadPos::hashLength + 1 - trimEnd; p++) {
    readPos.pos = p;
    valid = 1;
    for (p2 = 0; p2 < CountedReadPos::hashLength and valid; p2++) {
      if (unmasked_nuc_index[seqListPtr[curSeq].seq[p+p2]] >= 4)
				valid = 0;
    }
    if (valid and hashTable.Store(readPos, calcHashValue)) 
      newTupleStored = 1;
  }
  return newTupleStored;
}


void HashToReadPositionList(ReadPosHashTable &hashTable,
														SimpleSequenceList &sequences,
														ssize_t hashLength,
														ReadPositions &readPos ) {
	ssize_t i;
	ssize_t nReadPos = 0;
  ReadPosHashTable::iterator it, end;

  for (i = 0; i < hashTable.size; i++ ) {
    end = hashTable.table[i].End();
    it  = hashTable.table[i].Begin();
    while (it < end) {
      ++nReadPos;
      ++it;
    }
  }
	readPos.resize(nReadPos);
	ssize_t rp = 0;
	for (i = 0; i < hashTable.size; i++ ) {
    end = hashTable.table[i].End();
    it  = hashTable.table[i].Begin();
    while (it < end) {
			readPos[rp].read = (*it).read;
			readPos[rp].pos  = (*it).pos;
			rp++;
			++it;
		}
	}
	
	// Set up the functor to compare the tuples
  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &sequences;
  comp.length = hashLength;

	  // Do quicksort on the read positions
  std::sort(readPos.begin(), readPos.end(), comp);
}
