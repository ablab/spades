/***************************************************************************
 * Title:          MapReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "NumericHashedSpectrum.h"
#include "utils.h"
#include "IntegralTupleStatic.h"

using namespace std;

ssize_t hamming(char *c1, char *c2, ssize_t len) {
	char *cend;
	cend = c1 + len;
	ssize_t h = 0;
	for (; c1 != cend; c1++, c2++) {
		if (*c1 != *c2) ++h;
	}
	return h;
}
	
int main(int argc, char* argv[]) {


	string genomeFileName, readsFileName;
	if (argc < 3) {
		cout << "usage: mapReads genome reads [-maxHamming maxHamming (3)]" << endl;
		exit(0);
	}
	genomeFileName = argv[1];
	readsFileName  = argv[2];
	int argi = 3;
	ssize_t maxHammingDist = 3;
	while (argi < argc){ 
		if (strcmp(argv[argi], "-maxHamming") == 0) {
			maxHammingDist = atoi(argv[++argi]);
		}
		++argi;
	}

	DNASequence genome, genomeRC;
	ifstream readsIn;
	openck(readsFileName, readsIn, std::ios::in);

	SeqReader::GetSeq(genomeFileName, genome);
	
	MakeRC(genome, genomeRC);
	
	int vertexSize = 20;

	ssize_t dbSize = (genome.length - vertexSize + 1) * 2;
	CountedIntegralTuple *genomeDB = 
		new CountedIntegralTuple[dbSize];

	//	CountedIntegralTuple::tupleSize = vertexSize;
	CountedIntegralTuple::SetTupleSize(vertexSize);
	ssize_t i;
	for (i = 0; i < genome.length - vertexSize + 1; i++ ) {
		genomeDB[i].StringToTuple(&genome.seq[i]);
		genomeDB[i].count = i;
	}
	ssize_t irc = genome.length;
	for (i = 0; i < genome.length - vertexSize + 1; i++, irc++ ) {
		genomeDB[irc].StringToTuple(&genomeRC.seq[i]);
		genomeDB[irc].count = -i;
	}
	cerr << "sorting genomedb" << endl;
	std::sort(genomeDB,genomeDB + dbSize);
	cerr << "mapping reads" << endl;
	
	DNASequence read, readRC;
	ssize_t readFound;
	ssize_t readindex = 0;
	ssize_t numForSearched, numRevSearched;
	numForSearched = numRevSearched = 0;
	//UNUSED// ssize_t numRCSearched = 0;
	readindex = -1;
	while (SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		if (readindex % 10000 == 0) {
			cerr << readindex << endl;
		}

		++readindex;
		CountedIntegralTuple readTuple;
		readTuple.StringToTuple(read.seq);
		ssize_t readIndex = LookupBinaryTuple(genomeDB, dbSize, readTuple);
		readFound = 0;
		readRC.Reset(0);
		MakeRC(read, readRC);
		if (readIndex == -1) {
			++numRevSearched;

			readTuple.StringToTuple(readRC.seq);
			readIndex = LookupBinaryTuple(genomeDB, dbSize, readTuple);
			if (readIndex == -1)
				continue;

			// 
			// Make sure the match isn't in the forward strand.
			//
			if (genomeDB[readIndex].count < 0)
				continue;
			ssize_t h;
			while (readIndex < dbSize and 
						 genomeDB[readIndex] ==  readTuple // compare tuples
						 //						 genomeDB[readIndex].tuple ==  readTuple.tuple
						 ) {
				if (genomeDB[readIndex].count > 0 and
						(h = hamming((char*) &(genome.seq[genomeDB[readIndex].count]),
												 (char*) readRC.seq, read.length)) < maxHammingDist) {
					cout << genomeDB[readIndex].count << "\t" << h << " 1 " << readindex << endl;
				}
				++readIndex;
			}
		}
		else {
			//			cout << "looking up forward dir " << readIndex << endl;
			++numForSearched;
			if (genomeDB[readIndex].count < 0)
				continue;
			ssize_t h;

			while (readIndex < dbSize and
						 genomeDB[readIndex] == readTuple // compare tuples
						 //						 genomeDB[readIndex].tuple ==  readTuple.tuple
						 ) {
				
				if (genomeDB[readIndex].count >= 0 and
						(h = hamming((char*) &(genome.seq[genomeDB[readIndex].count]),
												 (char*) read.seq, read.length)) < maxHammingDist) {
					cout << genomeDB[readIndex].count << "\t" << h << " 0 " << readindex << endl;
				}
				else if (genomeDB[readIndex].count < 0 and
						(h = hamming((char*) &(genome.seq[genome.length + genomeDB[readIndex].count - read.length]),
												 (char*) readRC.seq, read.length)) < maxHammingDist) {
					cout << genome.length + genomeDB[readIndex].count - read.length + 1 
							 << "\t" << h << " 1 " << readindex << endl;
				}
				++readIndex;
			}
		}
	}
	cerr <<"Searched " << numForSearched << " forward and " << numRevSearched << " reverse." << endl;
	return 0;
}

