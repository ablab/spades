/***************************************************************************
 * Title:          BitSpectrum.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BIT_SPECTRUM_H_
#define BIT_SPECTRUM_H_

#include <string>
#include "assert.h"
#include <iostream>
#include <stdlib.h>
//#include <malloc.h>
#include "utils.h"
#include "Spectrum.h"

template<typename T>
class BitSpectrum : public Spectrum<T> {
 private:
	ssize_t GetMultiplicity(size_t pos, ssize_t slot) {/*
		std::cout << "gm: " << std::hex << (ssize_t) tupleList[pos] << " mo ";
		std::cout << slot << " " << std::oct << (ssize_t) includeMask[slot] << " slot  " << slot*2 << " maksed " ;
		std::cout << std::oct << (ssize_t) (tupleList[pos] & includeMask[slot]) ;
		std::cout << " " << (ssize_t) ((tupleList[pos] & includeMask[slot]) >> (slot*2)) << std::endl;
																					*/
		/*
		std::cout << (ssize_t) tupleList[pos] << " " << (ssize_t) includeMask[slot] ;
		std::cout << " " << (ssize_t) (tupleList[pos] & includeMask[slot]);
		std::cout << " " << (ssize_t) ((tupleList[pos] & includeMask[slot]) >> (slot * 2)) 
							<< std::endl;
		*/
		return (tupleList[pos] & includeMask[slot]) >> (slot * 2);
	}

 public:
	ssize_t findOnlySolid;
	ssize_t size() {
		return listLength;
	}
	BitSpectrum(int ts) {
		this->tupleSize = ts;
		//		T::tupleSize = ts;
		//		SetTupleSize(ts);
		T::SetTupleSize(ts);
		InitList();
		findOnlySolid = 0;
	}
	
	void FindOnlySolid() {
		findOnlySolid = 1;
	}

	BitSpectrum(std::string &spectrumFile) {
		Read(spectrumFile);
	}
	
	_SZT_ listLength;
	_SZT_ numTuples;
	unsigned char* tupleList;

	static unsigned char excludeMask[4];
	static unsigned char includeMask[4];

	void InitList() {
		numTuples = 1;
		numTuples <<= (this->tupleSize)*2;		

		//		std::cout << "numtuples: " << numTuples << std::endl;
		listLength = numTuples / 4;
		//		std::cout << "list length: " << listLength << std::endl;
		// Allocate space for the last tuple
		if (numTuples % 4 != 0) 
			listLength++;
		tupleList = new unsigned char[listLength];
	}


	ssize_t GetMultiplicity(size_t hashValue) {
		size_t pos = hashValue / 4;
		assert(pos < listLength);
		ssize_t slot = hashValue % 4;
		return GetMultiplicity(pos, slot);
	}

	void Write(std::string &fileName, ssize_t minMult=0) {
		// ignore the min mult for now.
		std::ofstream listOut;
		openck(fileName, listOut, std::ios::out);
		listOut << listLength << std::endl;
		listOut << numTuples  << std::endl;
		listOut.write((const char*) tupleList, listLength);
	}

	void Read(std::string &fileName, ssize_t minMult=0) {
		// ignore min mult for now
		std::ifstream listIn;
		openck(fileName, listIn, std::ios::in);
		listIn >> listLength;
		listIn >> numTuples;
		tupleList = new unsigned char[listLength];
		memset(tupleList, 0, listLength);
		std::string line;
		std::getline(listIn, line);
		listIn.read((char*) tupleList, listLength);
		//		if (!listIn.read((char*) tupleList, listLength)) {
		//			return 0;
		//		}
		//		else {
		//			return 1;
		//		}
	}

	ssize_t FindTuple(T &tuple) {
		size_t hashValue;
		if (tuple.GetHashValue(hashValue) < 0)
			return -1;
		ssize_t mult = GetMultiplicity(hashValue);
		if (findOnlySolid && mult >= 3)
			return 1;
		if (!findOnlySolid && mult != 0)
			return 1;
		
		// otherwise, the tuple isn't found.  It is 
		// either too low multiplicity and we are
		// only looking for solid, or there are no mult limits
		// and it is simply not found.
		return -1;
	}


	ssize_t IncrementMult(T &tuple) {
		size_t hashValue;
		ssize_t result = tuple.GetHashValue(hashValue);
		// Don't store hash values that have
		// masked nucleotides
		if (result < 0) return -1;
		size_t pos = hashValue / 4;
		ssize_t slot = hashValue % 4;
		assert(pos < listLength);
		ssize_t curMult = GetMultiplicity(pos, slot);
		//		std::cout << "location: "  << pos << " " << slot << " cur mult " << curMult << std::endl;
		if (curMult < 3) {
			curMult++;
			curMult <<= (slot*2);
			// clear the multiplicity at this spot
			tupleList[pos] &= excludeMask[slot];

			// put back in the new multiplicity
			tupleList[pos] |= curMult;
		}
		return curMult;
	}

	ssize_t CountSolid() {
		size_t pos;
		ssize_t slot;
		ssize_t numSolid = 0;
		for (pos = 0; pos < listLength-1; pos++) {
			for (slot = 0; slot < 4; slot++) {
				ssize_t mult = GetMultiplicity(pos, slot);
				if (mult >= 3) {
					numSolid++;
				}
			}
		}
		return numSolid;
	}

	ssize_t CountAll(size_t freq[]) {
		size_t pos, slot;
		ssize_t numSolid = 0;
		freq[0] = freq[1] = freq[2] = freq[3];
		for (pos = 0; pos < listLength-1; pos++) {
			for (slot = 0; slot < 4; slot++) {
				ssize_t mult = GetMultiplicity(pos, slot);
				freq[mult]++;
			}
		}
		ssize_t last = numTuples % 4;
		if (last == 0)
			last = 4;
		for (slot = 0; slot < last; slot++) {
			ssize_t mult = GetMultiplicity(pos,slot);
			if (mult >= 0) 
				freq[mult]++;
		}
		return numSolid;
	}

};


template<typename T>
unsigned char BitSpectrum<T>::excludeMask[4] = {0xfc, 0xf3, 0xcf, 0x3f};

template<typename T>
unsigned char BitSpectrum<T>::includeMask[4] = {0x03, 0x0c, 0x30, 0xc0};

#endif
