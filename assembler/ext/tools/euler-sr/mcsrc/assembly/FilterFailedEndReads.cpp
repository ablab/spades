/***************************************************************************
 * Title:          FilterFailedEndReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/27/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "IntegralTupleStatic.h"
#include "SeqUtils.h"

using namespace std;

void PrintUsage() {
	cout << "usage: filterFailedEndReads readsIn spectrumName tupleSize readsOut [-minMult M] [-minTuples t]" << endl;
}


int main(int argc, char* argv[]) {

	string readsFileName;
	string spectrumFileName;
	string readsOutName;
	int tupleSize;

	readsFileName = argv[1];
	spectrumFileName = argv[2];
	tupleSize = atoi(argv[3]);
	readsOutName  = argv[4];
	//	IntegralTuple::tupleSize = tupleSize;
	IntegralTuple::SetTupleSize(tupleSize);
	int argi = 5;
	ssize_t minTuples = 1;
	ssize_t minMult   = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minTuples") == 0) {
			minTuples = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minMult") == 0) {
			minMult = atoi(argv[++argi]);
		}
		else {
			cout << "Bad option: " << argv[argi] << endl;
			PrintUsage();
			exit(0);
		}
		++argi;
	}
	
	CountedIntegralTupleDict spectrum;
	spectrum.InitFromFile(spectrumFileName, minMult);

	ifstream readsIn;
	ofstream readsOut;
	openck(readsFileName, readsIn, std::ios::in);
	openck(readsOutName, readsOut, std::ios::out);

	DNASequence read, readRC;
	ssize_t numSkipped = 0;
	while (SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		MakeRC(read, readRC);

		// Try and find tuples in the RC in the directed tuple list.
		ssize_t revCount = 0;

		ssize_t i;
		CountedIntegralTuple tuple;
		// int index;

		for (i = 0; i < read.length - tupleSize + 1 and revCount < minTuples; i++) {
			if (tuple.StringToTuple(&(readRC.seq[i]))) {
				if (spectrum.DictLookupBinaryTuple(tuple) != -1) {
					++revCount;
				}
			}
		}
		if (revCount >= minTuples) {
			read.PrintlnSeq(readsOut);
		}
		else {
			/*
			if (read.length > tupleSize) {
				read.PrintlnSeq(cout);
			}
			*/
			++numSkipped;
		}
		read.Reset();
		readRC.Reset();
	}
	cout << "Filtered: " << numSkipped << " reads." << endl;
	return 0;
}
