/***************************************************************************
 * Title:          PrintUnmaskedContigs.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "utils.h"


int main(int argc, char* argv[]) {
	std::string inFileName, outFileName;
	inFileName = argv[1];
	outFileName = argv[2];
	
	std::ifstream inFile;
	std::ofstream outFile;
	openck(inFileName, inFile, std::ios::in);
	openck(outFileName, outFile, std::ios::out);
	DNASequence inSeq, contig;

	inSeq._masked = 1;
	SeqReader::MaskRepeats();
	while(SeqReader::GetSeq(inFile, inSeq, SeqReader::noConvert)) {
		ssize_t p;
		ssize_t contigStart = -1, contigEnd = -1;
		ssize_t maxGap = 30;
		ssize_t end;
		for (p = 1; p < inSeq.length; p++) {
			if (isupper(inSeq.seq[p]) and islower(inSeq.seq[p-1]))
				contigStart = p;
			if (islower(inSeq.seq[p]) and isupper(inSeq.seq[p-1])) {
				ssize_t p2;
				end= p + maxGap;
				if (end > inSeq.length) end = inSeq.length;
				ssize_t gapEnded = 0;
				for (p2 = p; p2 < end; p2++) {
					if (isupper(inSeq.seq[p2]))
						gapEnded = 1;
				}
				if (gapEnded)
					continue;

				contigEnd  = p;
			
				ssize_t length = contigEnd - contigStart;
				if (length > 1) {
					contig.seq = (unsigned char*) &inSeq.seq[contigStart];
					contig.length = length;
					std::stringstream namestrm;
					namestrm << contigStart << "_" << contigEnd;
					contig.namestr = namestrm.str();
					contig.PrintlnSeq(outFile);
				}
			}
		}
		// print the last contig.
		if (contigEnd != inSeq.length - 1) {
			contig.length = inSeq.length - contigStart;
			contig.seq = &inSeq.seq[contigStart];
			contig.PrintlnSeq(outFile);
		}
	}
	return 0;
}
	

