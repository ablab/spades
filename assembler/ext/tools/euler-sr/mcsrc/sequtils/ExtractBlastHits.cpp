/***************************************************************************
 * Title:          ExtractBlastHits.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "parseblast/BlastParser.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "ParseTitle.h"
#include "utils.h"
#include "compatibility.h"
#include <map>


int main(int argc, char* argv[]) {

	std::string seqFile, blastFile, readsFile, outSeqFile;

	if (argc != 5) {
		std::cout <<" usage: extractBlastHits genome blasttab reads readsOut" << std::endl;
		exit(1);
	}
	seqFile = argv[1];
	blastFile = argv[2];
	readsFile = argv[3];
	outSeqFile = argv[4];

	
	BlastResult blastResult;
  
  if (!(ReadBlastTable(blastFile, blastResult))) {
    std::cout << "Could not parse the blast file " << std::endl;
  }

	DNASequence ref;
	SeqReader::GetSeq(seqFile, ref);

	std::ifstream readIn;
	openck(readsFile, readIn, std::ios::in);
	DNASequence seq;
	std::map<std::string, ssize_t> readLengths;
	std::string fastaTitle;
	while (SeqReader::GetSeq(readIn, seq, SeqReader::noConvert)) {
		ParseTitle(seq.namestr, fastaTitle);
		readLengths[fastaTitle] = seq.length;
	}

	ssize_t br;
	_SSZT_ readLength; // TODO: see if could be size_t
	std::ofstream extractOut;
	DNASequence extractStr;
	extractStr._ascii = 1;
	openck(outSeqFile, extractOut, std::ios::out);
	ssize_t start, length;
	for (br = 0; br < blastResult.size(); br++ ) {
		if (blastResult[br]->hsps.size() > 0) {
			if (readLengths.find(blastResult[br]->queryName) != readLengths.end()) {
				readLength = readLengths[blastResult[br]->queryName];
				start = blastResult[br]->hsps[0].refPos -
					blastResult[br]->hsps[0].qryPos;
				
				length = std::min(readLength, ref.length - start);

				if (blastResult[br]->hsps[0].strand == 1) {
					start = ref.length - length - 1;
				}
				extractStr.seq = &ref.seq[start];
				extractStr.length = length;
				extractStr.namestr = blastResult[br]->queryName;
				extractStr.PrintSeq(extractOut);
				extractOut << std::endl;
			}
		}
	}
	return 0;
}
					

