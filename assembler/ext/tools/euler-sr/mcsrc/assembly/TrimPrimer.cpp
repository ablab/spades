/***************************************************************************
 * Title:          TrimPrimer.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "align/alignutils.h"
#include "IntegralTupleStatic.h"



int main(int argc, char* argv[]) {

	std::string seqFile, primerFile, seqOutFile;
	
	if (argc != 4) {
		std::cout << "usage: trimPrimer sequenceFile primerFile seqOutFile " << std::endl;
		return 0;
	}
	seqFile = argv[1];
	primerFile = argv[2];
	seqOutFile = argv[3];
	DNASequence primer, seq;
	
	SeqReader::GetSeq(primerFile, primer, SeqReader::noConvert);

	Score score(-1, 1, 5,1);
	
	//UNUSED// ssize_t s;
	ssize_t *locations;
	double alignScore;
	FloatMatrix scoreMat;
	//UNUSED// ssize_t l;
	DNASequence read;
	std::ifstream inFile;
	std::ofstream outFile;
	openck(seqFile, inFile, std::ios::in);
	openck(seqOutFile, outFile, std::ios::out);
	while (SeqReader::GetSeq(inFile, read, SeqReader::noConvert)) {
		locations  = NULL;
		alignScore = OverlapAlign(read, primer,  -1,1,3, locations, score.scoreMat);
		//		std::cout << read.namestr << " " << alignScore << std::endl;
		// may need to trim this sequence.
		ssize_t lastAlignPos = read.length - 1;

		if (locations[lastAlignPos]!= -1) {
			while(lastAlignPos > 0 and locations[lastAlignPos-1] != -1)
				--lastAlignPos;
		}
		ssize_t trimLength = read.length - lastAlignPos -1;
		//		std::cout << "trimming " << read.namestr << " by: " << trimLength << std::endl;
		read.length -= trimLength;
		read.PrintlnSeq(outFile);
		delete[] locations;
	}
	return 0;
}
