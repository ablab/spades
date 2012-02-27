/***************************************************************************
 * Title:          SplitFastaFile.cpp 
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


int main(int argc, char* argv[]) {

	std::string refFile;
	ssize_t maxSeqLength;
	std::string outFile;
	
	if (argc < 4) {
		std::cout << "usage: splitFastaFile refFile maxSeqLength outFile" << std::endl;
		exit(1);
	}
	refFile = argv[1];
	maxSeqLength = atoi(argv[2]);
	outFile = argv[3];

	std::ifstream in;
	openck(refFile, in);
	DNASequence seq, fragment;
	ssize_t pos;
	std::ofstream out;
	openck(outFile, out, std::ios::out);
	std::stringstream title;
	ssize_t index;
	while(SeqReader::GetSeq(in, seq, SeqReader::noConvert)) {
		pos = 0;
		index = 0;
		while (pos < seq.length) {
			fragment.seq = &seq.seq[pos];
			if (pos + maxSeqLength <= seq.length)
				fragment.length = maxSeqLength;
			else
				fragment.length = (seq.length - pos);
			title.str("");
			title << seq.namestr << "_" << index;
			fragment.namestr = title.str();
			fragment.PrintlnSeq(out);
			out << std::endl;
			pos += maxSeqLength;
			index++;
		}
	}			
	return 0;
}
