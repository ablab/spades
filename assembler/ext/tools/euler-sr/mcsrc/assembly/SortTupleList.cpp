/***************************************************************************
 * Title:          SortTupleList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "SimpleSequence.h"
#include "DeBruijnGraph.h"
#include <string>
#include <iostream>
#include <fstream>
#include "IntegralTupleStatic.h"


class CountedTuple : public SimpleSequence {
public:
  ssize_t count;
};

typedef std::vector<CountedTuple> CountedTupleList;
  
int main(int argc, char* argv[]) {
  std::string tupleListFile, sortedTupleListFile;
  if (argc < 3) {
    std::cout << "usage: sortTupleList tuplesIn tuplesOut [-ignoreCount] [-minMult m] [-lockFile file]" << std::endl;
    exit(1);
  }
  tupleListFile = argv[1];
  sortedTupleListFile = argv[2];
  ssize_t readCount = 1;
  int argi = 3;
  ssize_t minMult = 0;
	ssize_t doLock = 0;
	std::string lockFileName;
  while (argi < argc) {
    if (strcmp(argv[argi], "-ignoreCount") == 0) {
      readCount = 0;
    }
    else if (strcmp(argv[argi], "-minMult") == 0) {
      minMult = atoi(argv[++argi]);
    }
		else if (strcmp(argv[argi], "-lockFile") == 0) {
			doLock = 1;
			lockFileName = argv[++argi];
		}
    ++argi;
  }
  std::ifstream tupleListIn;
  std::ofstream tupleListOut;
  std::vector<CountedTuple> tuples;
  openck(tupleListFile, tupleListIn, std::ios::in);
	_INT_ lockFileDes;
  if (doLock){ 
		WaitLock(lockFileName, lockFileDes);
	}
  openck(sortedTupleListFile, tupleListOut, std::ios::out);
	
  ssize_t numTuples = 0;
  ssize_t count;
  if (minMult == 0) {
    // read it all
    if (! (tupleListIn >> numTuples)) {
      std::cout << "error, there should be a number of tuples at the" << std::endl;
      std::cout << "beginning of  " << tupleListFile << std::endl;
      exit(1);
    }
  }
  else {
    ssize_t numAboveThresh = 0;
    std::string st;
    ssize_t thresh;
    ssize_t oldNumTuples;
    tupleListIn >> oldNumTuples;
    while (tupleListIn) {
      if (!(tupleListIn >> st >> thresh))
				break;
      if (thresh >= minMult) {
				numAboveThresh++;
      }
    }
    numTuples = numAboveThresh;
    tupleListIn.close();
    tupleListIn.clear();
    openck(tupleListFile, tupleListIn, std::ios::in);
    // throw away the first line
    tupleListIn >> count;
  }

  tuples.resize(numTuples);
  std::vector<ReadPos> readPositions;
  readPositions.resize(numTuples);
  ssize_t i;
  std::string tuple;
  ssize_t pos = 0;
  while (tupleListIn) {
    if (!(tupleListIn >> tuple)) {
			break;
		}
    //std::cout << "got tuple: " << tuple << std::endl;
    if (readCount)
      tupleListIn >> count;
    
    if (!readCount or (count >= minMult)) {
      tuples[pos].seq = (unsigned char*) new unsigned char[tuple.size()];
      tuples[pos].length = tuple.size();
      memcpy((void*) tuples[pos].seq, (void*) tuple.c_str(), tuple.size());
      tuples[pos].count = count;
      
      readPositions[pos].read = pos;
      readPositions[pos].pos = 0;
      pos++;
    }
  }
  
	if (doLock) {
		ReleaseLock(lockFileDes);
	}

  CompareTuples<CountedTupleList> comp;
  comp.sequencesPtr = &tuples;
  comp.length = tuple.size();
  
  std::sort(readPositions.begin(), readPositions.end(), comp);
  unsigned char* tmpSeq = new unsigned char[tuple.size() + 1];
  tmpSeq[tuple.size()] = 0;

	if (doLock) {
		WaitLock(lockFileName, lockFileDes);
	}
  tupleListOut << readPositions.size() << std::endl;
  for (i = 0; i < readPositions.size(); i++ ) {
    memcpy(tmpSeq, tuples[readPositions[i].read].seq, tuple.size());
    tupleListOut << tmpSeq;
    if (readCount) 
      tupleListOut << " " << tuples[readPositions[i].read].count;
    tupleListOut << std::endl;
  }
  tupleListOut.close();
	
	if (doLock) {
		ReleaseLock(lockFileDes);
	}
  return 0;
}
