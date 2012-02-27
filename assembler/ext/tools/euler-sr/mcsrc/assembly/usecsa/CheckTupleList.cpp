/***************************************************************************
 * Title:          CheckTupleList.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/26/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "bbbwt/BBBWTQuery.h"
#include "../Spectrum.h"
#include "../Tuple.h"
#include "../ListSpectrum.h"
#include "../IntegralTupleStatic.h"

int main(int argc, char* argv[]) {
	std::string spectrumFile, csaFile;
	if (argc < 3) {
		std::cout <<"usage: checkTupleList spectrumFile csaFile [maxMult] " << std::endl;
		exit(0);
	}
	spectrumFile = argv[1];
	csaFile      = argv[2];
	ssize_t maxMult = 300;
	if (argc == 4) {
		maxMult = atoi(argv[3]);
	}
	std::vector<ssize_t> correctCount, incorrectCount;
	correctCount.resize(maxMult);
	incorrectCount.resize(maxMult);
	std::fill(correctCount.begin(), correctCount.end(), 0);
	std::fill(incorrectCount.begin(), incorrectCount.end(), 0);
	
	ListSpectrum<Tuple> spectrum;
	BBBWT csa;
	spectrum.Read(spectrumFile, 0);
	BW::Read(csaFile, csa);
	DNASequence query;
	ssize_t i;
	ssize_t nTuples = spectrum.size();
	ssize_t low, high;
	//UNUSED//	int tupleSize = spectrum.tupleSize;

	for (i = 0; i < nTuples; i++) {
		query.seq = (unsigned char*) spectrum[i].ToString();
		query.length = spectrum.tupleSize;
		BW::Query(query, csa, low, high);
		ssize_t mult = spectrum[i].GetMult();
		if (high - low > 0) {
			if (mult < maxMult-1) {
				correctCount[mult]++;
			}
			else {
				correctCount[maxMult-1]++;
			}
		}
		else {
			if (mult < maxMult-1) {
				incorrectCount[mult]++;
			}
			else {
				incorrectCount[maxMult-1]++;
			}
		}
	}
	for (i = 0; i < maxMult; i++) {
		std::cout << correctCount[i] << " " << incorrectCount[i] << std::endl;
	}
	return 0;
}
	
