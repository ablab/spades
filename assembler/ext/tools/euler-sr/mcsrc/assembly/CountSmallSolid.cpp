/***************************************************************************
 * Title:          CountSmallSolid.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BitSpectrum.h"
#include "StringTuple.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {

	std::string spectrumFile;
	if (argc < 2) {
		std::cout << "usage: countSmallSolid file" << std::endl;
		exit(0);
	}

	int argi = 1;
	spectrumFile = argv[argi];
	BitSpectrum<StringTuple> spectrum(spectrumFile);
	size_t freq[4];
	freq[0] = freq[1] = freq[2] = freq[3] = 0;
	std::cout << "total solid: " << spectrum.CountSolid() << std::endl;
	
	spectrum.CountAll(freq);
	std::cout << freq[0] << " " << freq[1] << " "
						<< freq[2] << " " << freq[3] << std::endl;
	return 0;
}
