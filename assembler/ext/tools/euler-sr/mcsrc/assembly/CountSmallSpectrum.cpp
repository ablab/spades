/***************************************************************************
 * Title:          CountSmallSpectrum.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "BitSpectrum.h"
#include "StringTuple.h"
#include "BufferedSeqReader.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


#include <string>
#include <vector>


#define MAX_TUPLE_SIZE 16

void PrintUsage() {
	std::cout << "usage: countSmallSpectrum spectrumFile [fastaFiles ...] -tupleSize T -spectrumFile file" << std::endl;
	exit(0);
}

int main(int argc, char* argv[]) {

	std::vector<std::string> readFiles;
	std::string spectrumFile;
	int argi = 1;
	int tupleSize = 0;
	while (argi < argc) {
		if (IsOption(argv[argi])) {
			if (IsOption(argv[argi], "-tupleSize")) {
				tupleSize = atoi(argv[++argi]);
			}
			else if (IsOption(argv[argi], "-spectrumFile")){ 
				spectrumFile  = argv[++argi];
			}
			else {
				std::cout << "unknown option " << argv[argi];
			}
		}
		else {
			if (argi < argc) 
				readFiles.push_back(argv[argi]);
		}
		++argi;
	}

	if (tupleSize == 0) {
		PrintUsage();
		std::cout << "you must specify a tuple size" << std::endl;
	}
	if (tupleSize > MAX_TUPLE_SIZE) {
		std::cout << "maximum tuple size is " << MAX_TUPLE_SIZE << std::endl;
	}

	ssize_t r;
	BufferedSeqReader<1000> seqReader;
	DNASequence read;
	StringTuple tuple;
	tuple.tupleSize = tupleSize;
	
	BitSpectrum<StringTuple> spectrum(tupleSize);
	ssize_t readNumber = 0;
	DNASequence readRC;
	ssize_t numHashed = 0;
	ssize_t numReads = 0;
	ssize_t numPositions = 0;
	for (r = 0; r < readFiles.size(); r++ ){ 
		seqReader.Init(readFiles[r]);
		std::ifstream in;
		openck(readFiles[r], in, std::ios::in);
		while(SeqReader::GetSeq(in, read, SeqReader::noConvert)) {
		//		while(seqReader.GetSeq(read)) {
			//			std::cout << read.namestr << std::endl;
			numReads+=2;
			if (readNumber % 1000 == 999) {
				std::cout << "." << std::flush;
			}
			if (readNumber % 50000 == 49999) {
				std::cout << " " << (readNumber+1) << std::endl;
			}
				
			ssize_t p;
			ssize_t nTuples = read.length - tupleSize + 1;
			size_t hashValue;
			for (p = 0; p < nTuples; p++) {
				tuple.assign((char*) &read.seq[p]);
				spectrum.IncrementMult(tuple);
				numPositions++;
				if (tuple.GetHashValue(hashValue) >= 0) {
					numHashed++;
					//					std::cout << hashValue << std::endl;
				}
			}

			MakeRC(read, readRC);
			for (p = 0; p < nTuples; p++) {
				tuple.assign((char*) &readRC.seq[p]);
				spectrum.IncrementMult(tuple);
				numPositions++;
				if (tuple.GetHashValue(hashValue) >= 0)  {
					numHashed++;
					//					std::cout << hashValue << std::endl;
				}
			}
			readNumber++;
		}
		seqReader.Reset();
	}
	std::cout << std::endl;
	std::cout << " hashed:  " << numHashed << " of " << numPositions << " of " << numReads << std::endl;
	spectrum.Write(spectrumFile);
}
