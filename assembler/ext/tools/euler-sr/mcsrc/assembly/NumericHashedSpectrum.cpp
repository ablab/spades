/***************************************************************************
 * Title:          HashedSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashedSpectrum.h"




void HashedSpectrum::HashToReadPositionList(SimpleSequenceList &sequences,
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


HashedSpectrum::HashedSpectrum(HashValueFunctor &hashFunctionP) {
	hashFunction = hashFunctionP;
	hashTable.Init(hashFunction);
}

	
