/***************************************************************************
 * Title:          IntegralTuple.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

// This source file is used in library libassemble, so use IntegralTuple.h
// instead of IntegralTupleStatic.h
#include "IntegralTuple.h"
#include "utils.h"

//LongTuple IntegralTuple::MASK_OFF_FIRST = 0;
#include <iomanip>

void IntegralTuple::ToString() const {
		std::string st;
		ToString(st);
		std::cout << st << std::endl;
	}



ssize_t ReadMultBoundedBinaryTupleList(std::string &fileName,
																	 ssize_t minMult,
																	 std::vector<CountedIntegralTuple>&  tupleList) {
	std::ifstream in;
	openck(fileName, in, std::ios::in | std::ios::binary);

	ssize_t tupleSize_SSZT;
	in.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	CountedIntegralTuple::SetTupleSize(tupleSize_SSZT);

	ssize_t numTuples;
	in.read((char*) &numTuples, sizeof(ssize_t));
	std::cout << "reading through " << numTuples << std::endl;
	CountedIntegralTuple tuple;
	ssize_t numFreqTuples = 0;
	while(in) {
		if (in.read((char*) &tuple, sizeof(CountedIntegralTuple))) {
			if (tuple.count >= minMult) {
				tupleList.push_back(tuple);
				numFreqTuples++;
			}
		}
	}
	in.clear();
	in.seekg(std::ios::beg);

	std::cout << "about to store " << numFreqTuples << " tuples." << std::endl;
	tupleList.resize(numFreqTuples);


	in.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));

	ssize_t i = 0;
	in.read((char*) &numTuples, sizeof(ssize_t));
	while(in) {
		if (in.read((char*) &tuple, sizeof(CountedIntegralTuple))) {
			if (tuple.count >= minMult) {
				//tupleList.push_back(tuple);
				//				std::string tupleStr;
				//				tuple.ToString(tupleStr);
				//				tupleList[i].tuple = tuple.tuple;
				tupleList[i].CopyTuple(tuple);
				tupleList[i].count = tuple.count;
				i++;
			}
		}
	}
	return numFreqTuples;
}
