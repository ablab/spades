/***************************************************************************
 * Title:          BinSpectToAscii.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "IntegralTupleStatic.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	std::string binListName, asciiListName;
	int tupleSize;
	if (argc < 4) {
		std::cout << "usage: binSpectToAscii binListName tupleSize asciiListName" << std::endl;
		std::cout << "   -fasta        Prints the tuples as a fasta file so they" << std::endl
							<< "                 be assembled." << std::endl;
		std::cout << "   -printCount   Prints the count of each tuple." << std::endl;
		exit(1);
	}
	binListName = argv[1];
	tupleSize   = atoi(argv[2]);
	asciiListName = argv[3];
	ssize_t printFasta = 0;
	ssize_t printCount = 0;
	int argi = 4;
	while (argi < argc) {
		if(strcmp(argv[argi], "-fasta") == 0) {
			printFasta = 1;
		}
		else if (strcmp(argv[argi], "-printCount") == 0) {
			printCount = 1;
		}
		argi++;
	}
	
	std::ifstream binIn;
	std::ofstream asciiOut;
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);

	openck(binListName, binIn, std::ios::in | std::ios::binary);
	openck(asciiListName, asciiOut, std::ios::out);

	// TODO: Added tupleSize to file format; should remove it as a command-line parameter
	ssize_t tupleSize_SSZT;
	binIn.read((char*) &tupleSize_SSZT, sizeof(tupleSize_SSZT));
	if (tupleSize != tupleSize_SSZT) {
		std::cerr << "Warning: command-line parameter tupleSize=" << tupleSize
							<< " disagrees with binary file tupleSize=" << tupleSize_SSZT << std::endl;
	}
	if (printCount == 0) {
		IntegralTuple::SetTupleSize(tupleSize_SSZT);
	} else {
		CountedIntegralTuple::SetTupleSize(tupleSize_SSZT);
	}

	ssize_t nTuples;
	binIn.read((char*) &nTuples, sizeof(ssize_t));
	ssize_t i;
	IntegralTuple tuple;
	CountedIntegralTuple cTuple;
	std::string tupleStr;
	ssize_t seqNumber = 0;
	if (!printFasta) {
		asciiOut << tupleSize << std::endl;
		asciiOut << nTuples << std::endl;
	}

	for (i = 0; i < nTuples ; i++ ) {
		if (printFasta) {
			asciiOut << ">" << seqNumber << std::endl;
		}
		if (printCount == 0) {
			binIn.read((char*) &tuple, sizeof(IntegralTuple));
			tuple.ToString(tupleStr);
			asciiOut << tupleStr << std::endl;	
		}
		else {
			binIn.read((char*) &cTuple, sizeof(CountedIntegralTuple));
			cTuple.ToString(tupleStr);
		 asciiOut << tupleStr << " " << cTuple.count << std::endl;
		}
		seqNumber++;
	}
	return 0;
}
