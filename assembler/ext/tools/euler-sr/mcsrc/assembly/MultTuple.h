/***************************************************************************
 * Title:          MultTuple.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MULT_TUPLE_H_
#define MULT_TUPLE_H_

#include "DNASequence.h"
#include "SeqUtils.h"
#include "Tuple.h"

class MultTuple : public Tuple {
	public:
	ssize_t mult;
	ssize_t Valid() {
		ssize_t i;
		for ( i = 0; i < this->size(); i++) {
			if (numeric_nuc_index[(unsigned char) (*this)[i]] >= 4)
				return 0;
		}
		return 1;
	}
	ssize_t GetMult() {
		return mult;
	}
 MultTuple(const char *s) : Tuple(s) {
		mult = 1;
	}
	MultTuple() {}
	void assign(const char *charPtr, ssize_t length=tupleSize) {
		((std::string*)this)->assign((char*) charPtr, tupleSize);
		this->mult = 1;
	}
	ssize_t ReadLine(std::istream &in, ssize_t minMult=0) {
		if ((in >> *((std::string*) this) >> mult)) {
			std::string line;
			std::getline(in, line);
			if (mult >= minMult)
				return 1;
		}
		return 0;
	}
	ssize_t IncrementMult() {
		mult++;
		return mult;
	}
	MultTuple &operator=(const std::string &str) {
		*((std::string*) this) = str;
		return *this;
	}
};



#endif
