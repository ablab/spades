/***************************************************************************
 * Title:          AlignPairedSequences.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "Alignment.h"
#include "AlignmentPrinter.h"
#include "align/alignutils.h"

int main(int argc, char* argv[]) {

	std::string fileAName, fileBName;
	fileAName = argv[1];
	fileBName = argv[2];

	std::ifstream fileA, fileB;

	openck(fileAName, fileA, std::ios::in);
	openck(fileBName, fileB, std::ios::in);

	DNASequence seqA, seqB;
	if (! (SeqReader::GetSeq(fileA, seqA, SeqReader::noConvert))) {
		std::cout << "Could not read from " << fileAName << std::endl;
		return 1;
	}

	
	if (! (SeqReader::GetSeq(fileB, seqB, SeqReader::noConvert))) {
		std::cout << "Could not read from " << fileBName << std::endl;
		return 1;
	}

	Score scoreMat(-1, 1, 1, 1);

	while (1) {
		if (seqA.namestr == seqB.namestr) {
			// The two sequences have the same name, align them
			ssize_t *alignment;
			alignment = new ssize_t[seqA.length+1];
			IntMatrix score;
			IntMatrix path;
			CreateMatrix(score, seqA.length+1, seqB.length+1);
			CreateMatrix(path,  seqA.length+1, seqB.length+1);
			
			//			double alignScore = 
				Align(seqA, seqB, -1, 1, 5, alignment, scoreMat.scoreMat);

			PrintAlignment(seqA, seqB, 0, 0, alignment, seqA.length, std::cout);
			delete []alignment;
			
			if (!(SeqReader::GetSeq(fileA, seqA, SeqReader::noConvert)))
				// no more A, done.
				return 0;

			if (!(SeqReader::GetSeq(fileB, seqB, SeqReader::noConvert)))
				// no more B, done.
				return 0;

		}
		else if (seqA.namestr < seqB.namestr) {
			if (!(SeqReader::GetSeq(fileA, seqA, SeqReader::noConvert)))
				// no more A, done.
				return 0;
		}
		else if (seqB.namestr < seqA.namestr){ 
			if (!(SeqReader::GetSeq(fileB, seqB, SeqReader::noConvert)))
				// no more B, done.
				return 0;
		}
	}
	return 0;
}
		
