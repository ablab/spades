/***************************************************************************
 * Title:          StringListSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ListSpectrum.h"

template<typename T>
ssize_t ListSpectrum<T>::ReadSpectrum(std::string &fileName, ssize_t minMult) {
	std::ifstream kmerFile;
	openck(fileName, kmerFile, std::ios::in);
	ssize_t mult;
	ssize_t numPastThresh = 0;
	ssize_t numTuples;
	int tupleSize;
	if (!kmerFile >> numTuples) {
		std::cout << "Error, the spectrum must begin with a number." << std::endl;
		exit(1);
	}
		
	if (!kmerFile >> tupleSize) {
		std::cout << "Error, the spectrum must contain a tuple size." << std::endl;
		exit(1);
	}
	T tempTuple;
	//	T::tupleSize = 20;
	T::SetTupleSize(20); // TODO: check
		
	while(kmerFile) {
		T::ReadLine(kmerFile, tempTuple);
		if (tempTuple.GetMult() >= minMult and tempTuple.Valid()) 
			numPastThresh++;
	}
	kmerFile.close();
	kmerFile.clear();
	openck(fileName, kmerFile, std::ios::in);
	kmerFile >> numTuples;
		
	Resize(numPastThresh);
	ssize_t index = 0;
	while(kmerFile) {
		T::ReadLine(kmerFile, tempTuple);
		if (tempTuple.GetMult() >= minMult and tempTuple.Valid()) {
			tupleList[index] = tempTuple;
			index++;
		}
	}
}

template<typename T>
void ListSpectrum<T>::WriteSpectrum(std::string &fileName, 
																		ssize_t minMult) {
	ssize_t numSpectra = 0;
	ssize_t i;
	for (i = 0 ; i < size(); i++ ) {
		if (tupleList[i].GetMult() >= minMult) numSpectra++;
	}
	std::ofstream spectOut;
	openck(fileName, spectOut, std::ios::out);
	spectOut << numSpectra << std::endl;
	for (i = 0;i < tupleList.size(); i++) {
		if (tupleList[i].GetMult() >= minMult) {
			spectOut << tupleList[i] << std::endl;
		}
	}
	spectOut.close();
}

template<typename T>
ssize_t ListSpectrum<T>::IncrementKmer(T &tuple) {
	ssize_t listIndex;
	listIndex = FindKmer(tuple);
	if (listIndex != -1) {
		return tupleList[listIndex].IncrementMult();
	}
	return -1;
}

template <typename T>
ssize_t ListSpectrum<T>::FindKmer(char *valueStr) {
	T tempTuple(valueStr);
	return FindKmer(tempTuple);
	tempTuple.CleanUp();
}

template <typename T>
ssize_t ListSpectrum<T>::FindKmer(Tuple &tuple) {

  // Do a simple binary search for value in spectrum

  ssize_t minIndex = 0;
  ssize_t maxIndex = tupleList.size();

  ssize_t index = (maxIndex + minIndex) / 2;

  if (tuple < tupleList[0] or 
      tuple > tupleList[tupleList.size()-1]) {
    return -1;
  }

  while (index < maxIndex and tupleList [index] != tuple) {
    if (maxIndex -1 == minIndex) 
      // minIndex is as close to the value is as possible, so it 
      // is not going to be found.
      return -1;

    if (tupleList[index] < tuple) {
      minIndex = index;
    }
    else if (tuple < tupleList[index]) {
      maxIndex = index;
    }
    else {
      std::cout << "case that shouldn't exist!!!!\n";
      exit(0);
    }
    index = (maxIndex + minIndex) / 2;
  }
  return index;
}
