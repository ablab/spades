/***************************************************************************
 * Title:          NumericTuple.h
 * Author:         Mark Chaisson
 * Created:        2010
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NUMERIC_TUPLE_H_
#define NUMERIC_TUPLE_H_

// TODO: Move constant MAX_TUPLE_SIZE to compatibility.h
// TODO: Integrate with IntegralTuple.h
#define MAX_TUPLE_SIZE 31
#include <iostream>
#include <string>
#include <sstream>

#include "SeqUtils.h"

#define NucsPerByte (CHAR_BIT / 2)

class NumericTuple  {
 public:
	static int tupleSize;
	static ssize_t numBytes;
	//	char tuple[MAX_TUPLE_SIZE/4 + 1];
	char tuple[(MAX_TUPLE_SIZE + (NucsPerByte-1)) / NucsPerByte];
	ssize_t GetMult() {
		return 0;
	}

	static void SetTupleSize(int ts) {
		tupleSize = ts;
		numBytes = (ts + (NucsPerByte-1)) / NucsPerByte;
	}
	
	NumericTuple() {
	}

	ssize_t Length() {
		return tupleSize;
	}

	ssize_t Valid() {
		// This may only store valid tuples.
		return 1;
	}
	
	ssize_t ReadLine(std::istream &in, ssize_t minMult =0) {
		std::string tuple;
		if (!(in >> tuple)) {
			std::cout << "Error reading tuple." << std::endl;
			exit(1);
		}
		if (tupleSize == -1) {
			// determine the tuple size from the word that was read.
			SetTupleSize(tuple.size());
		}
		StringToTuple(tuple);
		std::string line;
		std::getline(in, line);
		// Parse the multiplicity
		std::stringstream linestrm(line);
		ssize_t mult = 0;
		linestrm >> mult;
		if (mult >= minMult) 
			return 1;
		else
			return 0;
	}

	void StringToTuple(const std::string &s) {
		unsigned char nuc;
		ssize_t bytePos = 0;
		ssize_t byteIndex = 0;
		ssize_t p;
		ssize_t b;
		for (b = 0; b < numBytes; b++) {
			tuple[b] =0 ;
		}
		for (p = 0; p < s.size(); p++ ) {
			nuc = unmasked_nuc_index[(unsigned char) s[p]];
			nuc <<= (2*bytePos);
			tuple[byteIndex] |= nuc;
			bytePos++;
			if (p % 4 == 3) {
				byteIndex++;
			}
			bytePos %= 4;
		}
	}

	bool operator<(const NumericTuple &rhs) const {
		ssize_t bytePos;
		ssize_t byteIndex;
		ssize_t i;
		unsigned char byte, rhsByte;
		byte = tuple[0];
		rhsByte = rhs.tuple[0];
		unsigned char nuc, rhsNuc;
		unsigned char maskFirst = 0x3;
		byteIndex = 0;
		for (i = 0; i < tupleSize; i++) {
			nuc    = nuc_char[byte    & maskFirst];
			rhsNuc = nuc_char[rhsByte & maskFirst];
			byte >>=2;
			rhsByte >>=2;

			if (nuc < rhsNuc)
				return 1;
			if (rhsNuc < nuc)
				return 0;

			bytePos++;
			if (i % 4 == 3) {
				byteIndex++;
				byte    = tuple[byteIndex];
				rhsByte = rhs.tuple[byteIndex];
			}
			bytePos %= 4;
		}

		// The two words are equal
		return 0;
	}

	bool operator>(const NumericTuple &rhs) const {
		return rhs < *this;
	}
	
	bool operator==(const NumericTuple &rhs) const {
		ssize_t b;
		for (b = 0; b < numBytes; b++ ){
			if (rhs.tuple[b] != tuple[b])
				return 0;
		}
		return 1;
	}
	
	bool operator!=(const NumericTuple &rhs) const {
		return ! (*this == rhs);
	}
	
	NumericTuple& operator=(const NumericTuple &rhs) {
		if (this == &rhs)
			return *this;

		ssize_t b;
		for (b = 0; b < numBytes; b++) 
			this->tuple[b] = rhs.tuple[b];
		return *this;
	}

	NumericTuple &operator=(const std::string &tupleString) {
		this->StringToTuple(tupleString);
		return *this;
	}
	
	void ToString(std::string &tupleString) const {
		tupleString.reserve(tupleSize);
		ssize_t byteIndex, bytePos;
		ssize_t i;
		byteIndex = 0;
		unsigned char nuc = tuple[byteIndex];
		for (i = 0; i < tupleSize; i++) {
			// translate binary into nuc
			tupleString[i] = nuc_char[nuc & 0x3];
			// 
			bytePos ++;
			if (i % 4 == 3) {
				++byteIndex;
				nuc = tuple[byteIndex];
			}
		}
	}

	friend std::ostream &operator<<(std::ostream &out, const NumericTuple &tup) {
		std::string tupleString;
		tup.ToString(tupleString);
		out << tupleString << " 0" << std::endl;
		return out;
	}

	friend std::istream &operator>>(std::istream &in, NumericTuple &tup) {
		tup.ReadLine(in);
		return in;
	}
	
	ssize_t IncrementMult() {
		// No-op
		return 0;
	}
	
};


/*
	unsigned char NumericTuple::maskIn[4]  = {0x03, 0x0c, 0x30, 0xc0};
	unsigned char NumericTuple::maskOut[4] = {0xfc, 0xf3, 0xcf, 0x3f};
*/


#endif
