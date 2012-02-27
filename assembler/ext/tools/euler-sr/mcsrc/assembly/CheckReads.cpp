/***************************************************************************
 * Title:          CheckReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "ReadMap.h"
#include "ParseTitle.h"
#include "IntegralTupleStatic.h"

using namespace std;

void PrintUsage() {
		std::cout << "usage: checkReads readsFile genomeFile " 
							<< std::endl;
		std::cout << "-map mapFile  Reads are truncated at map." << std::endl;
		std::cout << "-printBad file Print reads with errors to 'file'" << std::endl;
		std::cout << "-statusVector file Print 1/0 for correct/incorrect to 'file'" << endl;
}
int main(int argc, char* argv[]) {

	std::string readsFileName, genomeFileName;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	readsFileName = argv[1];
	genomeFileName = argv[2];
	
	int argi = 3;
	std::string mapFileName;
	ssize_t useMap = 0;
	ssize_t printBad = 0;
	string errReadFileName;
	ssize_t printStatusVector = 0;
	string statusVectorFileName;
	while (argi < argc) {
		if (strcmp(argv[argi], "-map") == 0) {
			mapFileName = argv[++argi];
			useMap = 1;
		}
		else if (strcmp(argv[argi], "-printBad") == 0) {
			printBad = 1;
			errReadFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-statusVector") == 0) {
			printStatusVector = 1;
			statusVectorFileName = argv[++argi];
		}
		else {
			PrintUsage();
			cout << "bad option: " << argv[argi] << endl;
			exit(0);
		}
		++argi;
	}

	DNASequence genome, genomeRC;
	SeqReader::GetSeq(genomeFileName,genome, SeqReader::noConvert);

	ReadMapList readMap;
	if (useMap) {
		ReadReadMapFile(mapFileName, readMap);
	}

	ssize_t nErrors = 0;
	std::ifstream readsIn;
	openck(readsFileName, readsIn, std::ios::in);
	std::ofstream errOut;
	if (printBad)
		openck(errReadFileName, errOut, std::ios::out);
	DNASequence read;
	ssize_t nTotal = 0;
	//UNUSED// ssize_t readStart, readEnd;
	ssize_t dir, pos;
	ssize_t readIndex = 0;
	//UNUSED// ssize_t readLength;
	std::ofstream statusVectorOut;
	if (printStatusVector) {
		openck(statusVectorFileName, statusVectorOut, std::ios::out);
	}
	while(SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
		if (! ParseKeyword(read.namestr, "strand", dir)) {
			std::cout << "Reads must have 'strand' value in FASTA title."<< std::endl;
			exit(0);
		}
		if (! ParseKeyword(read.namestr, "pos", pos)) {
			std::cout << "Reads must have 'strand' value in FASTA title."<< std::endl;
			exit(0);
		}
		
		if (useMap) {
			pos += readMap[readIndex].start;
		}
		
		//		std::cout << read.namestr << " : " << dir << " " << pos << std::endl;
		ssize_t p;
		ssize_t isError = 0;
		if (dir == 0) {
			for (p = 0; p < read.length; p++ ) {
				if (genome.seq[pos + p] != read.seq[p]) {
					read.seq[p] = tolower(read.seq[p]);
					isError++;
				}
			}
		}
		else {
			for (p = 0; p < read.length; p++) {
				if (genome.seq[pos + read.length - p - 1] != comp_ascii[read.seq[p]]) {
					read.seq[p] = tolower(read.seq[p]);
					isError++;
				}
			}
		}
		nErrors += isError;
		if (printStatusVector) {
			if (isError) {
				statusVectorOut << "0" << endl;
			}
			else {
				statusVectorOut << "1" << endl;
			}
		}
		if (isError and printBad) {
			read.PrintlnSeq(errOut);
		}
		nTotal += read.length;
	}

	double errorRate = 0;
	if (nTotal > 0) { errorRate = double(nErrors) / nTotal; }

	std::cout << nErrors << " " << nTotal << " " << errorRate << std::endl;

}
