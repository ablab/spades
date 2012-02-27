/***************************************************************************
 * Title:          SortedTupleList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SortedTupleList.h"




void MakeSortedTupleList(SimpleSequenceList& refSeq,
												 int tupleSize,
												 std::vector<CountedReadPos> & refPositions) {
	ssize_t refLength = refSeq[0].length;
	refPositions.resize((refLength - tupleSize)*2);
	ssize_t i;
	for (i = 0; i < refLength - tupleSize + 1; i++) {
		refPositions[i].read = 0;
		refPositions[i].pos = i;
		refPositions[i+refLength-tupleSize+1].read = 1;
		refPositions[i+refLength-tupleSize+1].pos = i;
	}
	
  CompareTuples<SimpleSequenceList> comp;
  comp.sequencesPtr = &refSeq;
  comp.length = tupleSize;

  std::sort(refPositions.begin(), refPositions.end(), comp);
}

void RemoveDuplicatedTuples(SimpleSequenceList &refSeqList,
														int tupleSize,
														std::vector<CountedReadPos> & refPositions) {
	
	ssize_t cur = 0;
	ssize_t i;
	CountedReadPos::sequences = &refSeqList;
	CountedReadPos::hashLength = tupleSize;
	for (i = 0; i < refPositions.size(); ) {
		refPositions[cur] = refPositions[i];
		i++;
		while(i < refPositions.size() and
					refPositions[cur] == refPositions[i]) {
			i++;
		}
		cur++;
	}

	refPositions.resize(cur);
}

		
ssize_t CompareTupleSeq(ReadPos &a, char* seq, SimpleSequenceList &sequences,
										int tupleSize, _INT_ print) {
  if (a.pos + tupleSize > sequences[a.read].length ) {
    std::cout << "at: " << a.read << std::endl;
    std::cout << "error, trying to access position: " << a.pos 
              << " should be no less than " << a.pos + tupleSize 
              << " < " << sequences[a.read].length << std::endl;
  }

  assert(a.pos + tupleSize <= sequences[a.read].length);
  // temp stuff
  if (print) {
    std::cout << "comparing str ";
    PrintTuple(sequences, a, tupleSize);
    std::cout << " and : ";
    char tmp[1000];
    strncpy(tmp, seq, tupleSize);
    tmp[tupleSize] = 0;
    std::cout << tmp << std::endl;
  }
  char *strptr;
  strptr = (char*) &(sequences[a.read].seq[a.pos]);
  ssize_t result = strncmp(strptr, seq, tupleSize);
  return result;
}

ssize_t CompareTupleSeq(ReadPos &a, ReadPos &b, 
										SimpleSequenceList &sequences,
										int tupleSize, _INT_ print) {
  assert(a.pos + tupleSize-1 < sequences[a.read].length);
  assert(b.pos + tupleSize-1 < sequences[b.read].length);
  if (print) {
    std::cout << "comparing: ";
    PrintTuple(sequences, a, tupleSize);
    std::cout << " and: ";
    PrintTuple(sequences, b, tupleSize);
    std::cout << std::endl;
  }
  return (strncmp((char*) &(sequences[a.read].seq[a.pos]),
		  (char*) &(sequences[b.read].seq[b.pos]),
		  tupleSize));
}

														
