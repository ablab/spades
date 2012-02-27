/***************************************************************************
 * Title:          TrimReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/10/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"

void PrintUsage() {
	std::cout << "trimReads in out [-trimFront t] [-trimEnd t]" << std::endl;
	exit(0);
}
int main(int argc, char* argv[]) {

	std::string inFile, outFile;
	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	inFile = argv[1];
	outFile = argv[2];
	std::ifstream in;
	std::ofstream out;
	openck(inFile, in, std::ios::in);
	openck(outFile, out, std::ios::out);
	int argi = 3;
	ssize_t trimFront = 0;
	ssize_t trimEnd   = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-trimFront") == 0) {
			trimFront = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-trimEnd") == 0) {
			trimEnd = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			exit(0);
		}
		++argi;
	}
		
	DNASequence seq;
	while (SeqReader::GetSeq(in, seq, SeqReader::noConvert)) {
		seq.length -= (trimFront + trimEnd);
		if (seq.length > 0) {
			seq.seq = &seq.seq[trimFront];
			seq.PrintlnSeq(out);
			seq.seq = &seq.seq[-trimFront];
		}
	}
	return 0;
}
	
