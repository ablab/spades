/***************************************************************************
 * Title:          NumericMultTuple.h
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NUMERIC_MULT_TUPLE_H_
#define NUMERIC_MULT_TUPLE_H_


//UNUSED// #define MAX_TUPLE_SIZE 30
#include <iostream>
#include <string>
#include <sstream>

#include "SeqUtils.h"
#include "NumericTuple.h"

class NumericMultTuple : public NumericTuple {
 public:
	ssize_t mult;
	ssize_t GetMult() {
		return mult;
	}
	
 NumericMultTuple() : NumericTuple() {
		mult = 0;
	}
	
	
	ssize_t ReadLine(std::istream &in, ssize_t minMult =0) {
		std::string tuple;
		if (!(in >> tuple >> mult)) {
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

		// Check multiplicity
		if (mult >= minMult) 
			return 1;
		else
			return 0;
	}


	
	ssize_t IncrementMult() {
		// No-op
		mult++;
		return mult;
	}
	
};

#endif
