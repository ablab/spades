/***************************************************************************
 * Title:          SortSequencesByName.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"


class CompareDNAPtrsByName {
public:
	ssize_t operator()(DNASequence *seqa, DNASequence *seqb) {
		return seqa->namestr < seqb->namestr;
	}
};


int main(int argc, char* argv[]) {
	
	if (argc < 3) {
		std::cout << "Sort sequences by fasta title." << std::endl;
		std::cout << "usage: sortSequences inFile outFile"<<std::endl;
		return 1;
	}
	std::string inFileName, outFileName;
	inFileName = argv[1];
	outFileName = argv[2];
	
	DNASequenceList inReads;
	ReadDNASequences(inFileName, inReads);
	std::vector<DNASequence*> readPtrList;
	readPtrList.resize(inReads.size());
	ssize_t i;
	for (i = 0; i < inReads.size(); i++) {
		readPtrList[i] = (DNASequence*) (&inReads[0] + i);
	}

	std::sort(readPtrList.begin(), readPtrList.end(), CompareDNAPtrsByName());

	std::ofstream outFile;
	openck(outFileName, outFile, std::ios::out);
	for (i = 0; i < inReads.size(); i++) {
		readPtrList[i]->PrintlnSeq(outFile);
		outFile << std::endl;
	}
	return 0;
}
