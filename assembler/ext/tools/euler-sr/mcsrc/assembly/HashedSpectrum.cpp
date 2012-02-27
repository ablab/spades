/***************************************************************************
 * Title:          HashedSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/18/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "HashedSpectrum.h"
#include "utils.h"





ssize_t HashedSpectrum::HashSequence(SimpleSequenceList &seqListPtr, ssize_t curSeq,
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
    if (valid and hashTable.Store(readPos, hashFunction)) 
      newTupleStored = 1;
  }
  return newTupleStored;
}

void HashedSpectrum::StoreSpectrum(SimpleSequenceList &seqList,
																	 ssize_t trimFront, ssize_t trimEnd ) {
	ssize_t s;
	ssize_t storeSeq;
	ssize_t seqNumber = 0;

	for (s = 0; s < seqList.size(); s++ ) {
		storeSeq = HashSequence(seqList, s, trimFront, trimEnd);
			
		seqNumber++;
		PrintStatus(seqNumber, 100000);
	}
}

void HashedSpectrum::HashToReadPositionList(SimpleSequenceList &sequences,
																						ReadPositions &readPos ) {
	ssize_t i;
	ssize_t nReadPos = 0;
  ReadPosHashTable::iterator it, end;
	ReadPosHashTable::HashPage *page;
	ReadPosHashTable::Data *data, *pageEnd;
  for (i = 0; i < hashTable.size; i++ ) {

		page = hashTable.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {
				++nReadPos;
			}
			page = page->next;
		}
		/*    end = hashTable.table[i].End();
    it  = hashTable.table[i].Begin();
    while (it < end) {
      ++nReadPos;
      ++it;
    }
		*/
  }
	readPos.resize(nReadPos);
	ssize_t rp = 0;
	for (i = 0; i < hashTable.size; i++ ) {
		page = hashTable.table[i].head;
		while (page != NULL) {
			pageEnd = &(*page).page[page->size];
			for (data = &(*page).page[0]; data != pageEnd; ++data) {

				/*    end = hashTable.table[i].End();
    it  = hashTable.table[i].Begin();
    while (it < end) {
				*/
				readPos[rp].read = data->read; //(*it).read;
				readPos[rp].pos  = data->pos;  // (*it).pos;
				rp++;
				++it;
			}
			page = page->next;
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

	
