/***************************************************************************
 * Title:          ShortSpectrumToLong.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "BitSpectrum.h"
#include "StringTuple.h"
#include "IntegralTupleStatic.h"

#include <iostream>


int main(int argc, char* argv[]) {

	std::string bitSpectrumFileName = "";
	bitSpectrumFileName = argv[1];
	int tupleSize = atoi(argv[2]);
	BitSpectrum<StringTuple> spectrum(tupleSize);
		
	spectrum.Read(bitSpectrumFileName);

	_SZT_ i;
	size_t t, tup;
	char* tuple = new char[tupleSize+1];
	tuple[tupleSize] = '\0';
	size_t mask = 0x3;
	ssize_t mult;
	size_t freq[4];
	spectrum.CountAll(freq);
	ssize_t nonZero = freq[1] + freq[2] + freq[3];
	std::cout << nonZero << std::endl;
	ssize_t numPrinted = 0;
	for (i = 0; i < spectrum.numTuples; i++) {
		tup = (size_t) i;
		for (t = 0; t < tupleSize; t++) {
			tuple[tupleSize - t - 1] = nuc_char[tup & mask];
			tup >>=2;
		}
		if ((mult = spectrum.GetMultiplicity(i)) > 0) {
			std::cout << tuple << " " << mult << std::endl;
			//			numPrinted++;
		}
		for (t = 0; t < tupleSize; t++) 
			tuple[t] = ' ';
	}
	std::cout << "would have printed" << numPrinted << std::endl;
}
