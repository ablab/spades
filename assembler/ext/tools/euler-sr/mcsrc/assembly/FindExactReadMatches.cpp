/***************************************************************************
 * Title:          FindExactReadMatches.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Tuple.h"
#include "ListSpectrum.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {
	
	std::string genomeSpectrumFileName, readsFileName;
	std::string containedFileName, notContainedFileName;
	if (argc < 5) {
		std::cout << "usage: findExactReadMatches reads genomeSpectrum contained notContained" 
							<< std::endl;
		exit(0);
	}
	readsFileName          = argv[1];
	genomeSpectrumFileName = argv[2];
	containedFileName      = argv[3];
	notContainedFileName   = argv[4];

	ListSpectrum<Tuple> listSpectrum;

	listSpectrum.Read(genomeSpectrumFileName, 0);

	DNASequence read;
	std::ifstream readIn;
	openck(readsFileName, readIn, std::ios::in);
	std::ofstream contOut, notContOut;
	openck(containedFileName, contOut, std::ios::out);
	openck(notContainedFileName, notContOut, std::ios::out);
	while(SeqReader::GetSeq(readIn, read, SeqReader::noConvert)) {
		//		Tuple::tupleSize = read.length;		
		Tuple::SetTupleSize(read.length);
		Tuple tempTuple;
		tempTuple.assign((char*)read.seq);
		ssize_t result;
		result = listSpectrum.FindTuple(tempTuple);
		//		std::cout << result << std::endl;
		if ( result == -1) {
			read.PrintSeq(notContOut);
			notContOut << std::endl;
		}
		else {
			read.PrintSeq(contOut);
			contOut << std::endl;
		}
	}
	return 0;
}
