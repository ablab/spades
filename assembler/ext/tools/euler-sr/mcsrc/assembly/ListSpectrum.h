/***************************************************************************
 * Title:          ListSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef LIST_SPECTRUM_H_
#define LIST_SPECTRUM_H_

#include "Spectrum.h"
#include "utils.h"

// TODO: Fix sparc compiler warning
//  "warning: `class ListSpectrum<T>' has virtual functions but non-virtual destructor"
// (T ranges over MultTuple, NumericTuple, StringMultTuple, StringTuple, Tuple)
template <typename T>
class ListSpectrum : public Spectrum<T> {

 public:
	std::vector<T> tupleList;

	ssize_t size() {
		return tupleList.size();
	}

	void Resize(ssize_t size) {
		tupleList.resize(size);
	}

	T &operator[](const _SSZT_ i) {
		return tupleList[i];
	}

	void Read(std::string &fileName, ssize_t minMult = 0) {
		std::ifstream kmerFile;
		openck(fileName, kmerFile, std::ios::in);
		//UNUSED// ssize_t mult;
		ssize_t numPastThresh = 0;
		ssize_t numTuples;
		if (!(kmerFile >> numTuples)) {
			std::cout << "Error, the spectrum must begin with a number." << std::endl;
			exit(1);
		}
		
		T tempTuple;
		ssize_t line = 0;
		while(kmerFile && line < numTuples) {
			//			tempTuple.ReadLine(kmerFile);
			line++;
			if (tempTuple.ReadLine(kmerFile, minMult)) {
				if (tempTuple.Valid()) {
					numPastThresh++;
				}
			}
		}
		kmerFile.close();
		kmerFile.clear();
		openck(fileName, kmerFile, std::ios::in);
		kmerFile >> numTuples;
		Resize(numPastThresh);
		ssize_t index = 0;
		line = 0;
		while(kmerFile and line < numTuples) {
			if (tempTuple.ReadLine(kmerFile, minMult)) {
				if (tempTuple.Valid()) {
					//					std::cout << "storing tuple at " << index << std::endl;
					//					Pause();
					tupleList[index] = tempTuple;
					index++;
				}
			}
			line++;
		}
		this->tupleSize = tempTuple.tupleSize;
		index++;
	}

	void Write(std::string &fileName, ssize_t minMult = 0) {
		ssize_t numSpectra = 0;
		ssize_t i;
		for (i =0 ; i < tupleList.size(); i++ ) {
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

	ssize_t IncrementMult(T &tuple) {
		ssize_t listIndex;
		listIndex = FindTuple(tuple);
		if (listIndex != -1) {
			return tupleList[listIndex].IncrementMult();
		}
		return 0;
	}

	ssize_t FindTuple(std::string &tupleStr) {
		T tempTuple;
		tempTuple = tupleStr;
		return FindTuple(tempTuple);

	}
	ssize_t FindTuple(T &tuple) {

		// Do a simple binary search for value in spectrum
		ssize_t minIndex = 0;
		ssize_t maxIndex = tupleList.size();

		ssize_t index = (maxIndex + minIndex) / 2;

		if (tuple < tupleList[0] or 
				tuple > tupleList[tupleList.size()-1]) {
			return -1;
		}

		while (index < maxIndex and tupleList[index] != tuple ) {
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

};

#endif
