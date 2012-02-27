/***************************************************************************
 * Title:          SpliceReads.cpp 
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


void PrintUsage() {
	std::cout << "usage: spliceReads prefixFile suffixFile outFile " << std::endl;
	std::cout << "   [-prefix prefixLength] " << std::endl
						<< "   [-suffix suffixLength] " << std::endl;
}

int main(int argc, char* argv[]) {
	std::string prefixReadsFile, suffixReadsFile;
	ssize_t prefixLength, suffixLength;
	
	std::string outFileName;
	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	prefixReadsFile = argv[argi++];
	suffixReadsFile = argv[argi++];
	outFileName     = argv[argi++];
	
	prefixLength = -1;
	suffixLength = -1;
	while (argi < argc) {
		if (strcmp(argv[argi], "-prefix") == 0) {
			prefixLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-suffix") == 0) {
			suffixLength = atoi(argv[++argi]);
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}

	std::ifstream prefixIn, suffixIn;
	std::ofstream readsOut;
	
	openck(prefixReadsFile, prefixIn, std::ios::in);
	openck(suffixReadsFile, suffixIn, std::ios::in);
	openck(outFileName, readsOut, std::ios::out);

	DNASequence prefix, suffix;
	DNASequence splice;
	while(SeqReader::GetSeq(prefixIn, prefix, SeqReader::noConvert)) {
		if (!SeqReader::GetSeq(suffixIn, suffix, SeqReader::noConvert)) {
			std::cout << "Error, the prefix file and suffix file must have "
								<< "the same number of reads." << std::endl;
			exit(1);
		}
		if (prefix.namestr != suffix.namestr) {
			std::cout << "ERROR! the prefix and suffix reads must have the " 
								<< "same FASTA name." << std::endl;
			exit(0);
		}
		ssize_t pre, suf;
		if (prefixLength == -1) 
			pre = prefix.length;
		else 
			pre = prefixLength;
		
		if (suffixLength == -1)
			suf = suffix.length;
		else
			suf = suffixLength;

		splice.Reset(pre + suf);
		
		// make sure there are not bounds overwrites
		assert(pre <= prefix.length);
		assert(suf <= suffix.length);

		memcpy(splice.seq, prefix.seq, pre);
		memcpy(splice.seq + pre, suffix.seq + suffix.length - suf, suf);
		
		splice.namestr= prefix.namestr;
		splice.PrintlnSeq(readsOut);
	}
	readsOut.close();
	return 0;
}

		
		
