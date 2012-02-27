/***************************************************************************
 * Title:          SortedTupleList.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SORTED_TUPLE_LIST_H_
#define SORTED_TUPLE_LIST_H_

#include "DNASequence.h"
#include "SimpleSequence.h"
#include "Spectrum.h"
#include "ReadPos.h"
#include "compatibility.h"


void MakeSortedTupleList(SimpleSequenceList& refSeq,
												 int tupleSize,
												 std::vector<CountedReadPos> & refPositions);

void RemoveDuplicatedTuples(SimpleSequenceList &refSeq,
														int tupleSize,
														std::vector<CountedReadPos> & refPositions);

ssize_t CompareTupleSeq(ReadPos &a, ReadPos &b, SimpleSequenceList &sequences,
										int tupleSize, _INT_ print = 0);

ssize_t CompareTupleSeq(ReadPos &a, char* seq, SimpleSequenceList &sequences,
										int tupleSize, _INT_ print = 0);


template<class T>
ssize_t LocateFirstTuple(SimpleSequenceList &sequences,
										 std::vector<T> &tuples,
										 int tupleSize,
										 char * tuple, _INT_ print = 0) {
  ssize_t index;
  index = LocateTuple(sequences, tuples, tupleSize, tuple);
  char* tuplePtr;
	if (index == -1)
		return -1;
  tuplePtr = (char*) &(sequences[tuples[index].read].seq[tuples[index].pos]);
  while (index >= 1 and 
				 CompareTupleSeq(tuples[index-1], tuplePtr, sequences, tupleSize)==0) 
    index--;
  return index;
}

template<class T>
ssize_t LocateTuple(SimpleSequenceList &sequences,
								std::vector<T> &tuples,
								int tupleSize,
								char * tuple, _INT_ print = 0) {
  ssize_t beg, end, cur;
  beg = 0;
  end = tuples.size();
  cur = (end + beg) / 2;
  
  //  std::cout << "size of list: "<< tuples.size() << std::endl;
  ssize_t comp = -1; // init to nonzero, so if while() doesn't execute then we return -1.
  ssize_t prev = -1;
  while ( beg < end  and 
					( comp = CompareTupleSeq(tuples[cur], tuple, sequences, tupleSize)) != 0 and
					prev != cur) {
    prev = cur;
    
    if (comp > 0) {
      end = cur;
    }
    else if (comp < 0) {
      beg = cur;
    }
    cur = (end + beg) / 2;
  }
  if (comp == 0) {
    return cur;
  }
  else {
    return -1;
  }
}


#endif
