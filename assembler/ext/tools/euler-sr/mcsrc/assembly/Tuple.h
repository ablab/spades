/***************************************************************************
 * Title:          Tuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef TUPLE_H_
#define TUPLE_H_

#include "DNASequence.h"
#include "SeqUtils.h"

#include <sstream>
class Tuple : public std::string {
 public:
	static int tupleSize;
	ssize_t Valid() {
		ssize_t i;
		for ( i = 0; i < this->size(); i++) {
			if (numeric_nuc_index[(unsigned char) (*this)[i]] >= 4)
				return 0;
		}
		return 1;
	}
	Tuple(const char* st) {
		std::string(s);
	}
	Tuple() {}
	ssize_t GetMult() {
		return 1;
	}
	static void SetTupleSize(int ts) {
		tupleSize = ts;
	}
	std::istream &ReadLine(std::istream &in, ssize_t minMult = 0) {
		in >> (*this);
		std::string line;
		std::getline(in, line);
		std::stringstream linestrm(line);
		if (tupleSize == 0) {
			tupleSize = this->size();
		}
		ssize_t mult;
		linestrm >> mult;
		return in;
	}
	const char *ToString() {
		return (const char*) this->c_str();
	}
	ssize_t IncrementMult() {
		return 0;
	}
	Tuple &operator=(const std::string &t) {
		*((std::string*)this) = t;
		return *this;
	}
	ssize_t GetHashValue(size_t &hashValue) {
		ssize_t i;
		if (this->size() == 0) {
			return -1;
		}
		hashValue = 0;
		// method 1
		for (i = 0; i < tupleSize; i++) {
			hashValue <<=2;
			if (unmasked_nuc_index[(unsigned char) (*this)[i]] >= 4)
				return -1;
			else
				hashValue += unmasked_nuc_index[(unsigned char) (*this)[i]];
		}
		return 1;
	}
};

template <typename T>
void Concatenate( T pre, unsigned char nuc, T& result) {
  // concatenate nuc to pre
  ssize_t length = pre.size();
	std::string resultStr;
  resultStr = "";
  resultStr.insert(0,(const char*) pre.ToString() + 1,length-1);
  resultStr.push_back(nuc_char[nuc]);
	result.assign((char*) resultStr.c_str());
  //  std::cout << "concatendated " << pre << " " << nuc_char[nuc] << "  to make: " << result << std::endl;
}

void InsertTuple(char* tuple, ssize_t tupleLen, ssize_t pos, unsigned char nuc);
void DeleteTuple(char* tuple, ssize_t tupleLen, ssize_t pos, unsigned char fill);
void MutateTuple(char* tuple, ssize_t pos, unsigned char nuc);


#endif
