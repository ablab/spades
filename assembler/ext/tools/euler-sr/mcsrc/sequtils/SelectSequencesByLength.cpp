/***************************************************************************
 * Title:          SelectSequencesByLength.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"

int main(int argc, char *argv[]) {
	std::string inFile, outFile;
	ssize_t minLength,  maxLength;

	if (argc < 3) {
		std::cout << "usage: selectSeq inFile outFile -s minLength (0) -m maxLength (infty) " << std::endl;
		exit(0);
	}

	inFile = argv[1];
	outFile = argv[2];

	int argi = 3;
	minLength = -1;
	maxLength = -1;
	while (argi < argc) {
		if (strcmp(argv[argi], "-s")==0) {
			minLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-m")==0) {
			maxLength = atoi(argv[++argi]);
		}
		++argi;
	}
	std::ifstream in;
	std::ofstream out;

	openck(inFile, in, std::ios::in);
	openck(outFile, out, std::ios::out);

	DNASequence seq;
	while (SeqReader::GetSeq(in, seq)) {
		if ((minLength == -1 or seq.length > minLength) and
				(maxLength == -1 or seq.length < maxLength)) {
			seq.PrintSeq(out);
			out << std::endl;
		}
	}
	return 0;
}
