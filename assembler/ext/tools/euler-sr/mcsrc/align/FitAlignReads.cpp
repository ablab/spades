/***************************************************************************
 * Title:          FitAlignReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "mctypes.h"
#include "utils.h"
#include "math.h"
#include "BufferedSeqReader.h"
#include "Alignment.h"
#include "AlignmentPrinter.h"
#include "align/alignutils.h"

void PrintUsage() {
		std::cout << "usage: fitAlignReads refFile readsFile " << std::endl;
		std::cout << " [-biasStart s] (1) bias mismatch on low-quality regions by s" << std::endl;
		std::cout << " [-biasEnd   s] (3) bias mismatch on high-quality regions by s"<< std::endl;
		std::cout << " [-gapped]          Perform gapped alignments" << std::endl;
		std::cout << " [-errProf  file] (none) Output error profile to 'file'" << std::endl;
}

int main(int argc, char* argv[]) {


	std::string refSeqFile, readsFile;

	if (argc < 3) {
		PrintUsage();
		exit(0);
	}
	refSeqFile = argv[1];
	readsFile  = argv[2];

	int argi = 3;
	ssize_t biasStart = 1;
	ssize_t biasEnd   = 3;
	ssize_t gappedAlign = 0;
	ssize_t printErrorProf = 0;
	std::string errProfName;
	while (argi < argc) {
		if (strcmp(argv[argi], "-biasStart") ==0) {
			biasStart = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-biasEnd") ==0) {
			biasEnd   = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-gapped" ) == 0) {
			gappedAlign = 1;
		}
		else if (strcmp(argv[argi], "-errProf" ) == 0) {
			errProfName = argv[++argi];
			printErrorProf = 1;
		}
		else {
			PrintUsage();
			std::cout << "bad option:" << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}
	
	DNASequence ref;
	ref._ascii = 1;
	SeqReader::GetSeq(refSeqFile, ref, SeqReader::noConvert);
 // BufferedSeqReader<1000> seqReader;
  // seqReader.Init(readsFile);
	std::ifstream readsIn;
	openck(readsFile, readsIn);
	DNASequence read;
	Score scoreMat(-1, 1, 1, 1);
	if (!SeqReader::GetSeq(readsIn, read, SeqReader::noConvert)) {
 // if (!seqReader.GetSeq(read)) {
		return 0;
	}

	std::ofstream errProfOut;
	
	if (printErrorProf) {
		openck(errProfName, errProfOut, std::ios::out);
	}

	
	ssize_t nRows = read.length + 1;
	//UNUSED// ssize_t nCols = ref.length;

	IntMatrix matchMat;
	CreateMatrix(matchMat, 4, 4);

	matchMat[0][0] = 1; matchMat[1][1] = 1; matchMat[2][2] = 1; matchMat[3][3] = 1;

	matchMat[0][1] = -1; matchMat[0][2] = -1; matchMat[0][3] = -1;
	matchMat[1][0] = -1; matchMat[1][2] = -1;	matchMat[1][3] = -1;
	matchMat[2][0] = -1; matchMat[2][1] = -1; matchMat[2][3] = -1;
	matchMat[3][0] = -1; matchMat[3][1] = -1;	matchMat[3][2] = -1;

	ssize_t readLength = read.length;
	ssize_t refLength  = ref.length;
	
	std::vector<double> bias;
	ssize_t p;
	bias.resize(readLength);
	for (p = 0; p < readLength; p++) {
		bias[p] = biasStart + (biasEnd - biasStart) * double(readLength - p)/nRows;
		//		bias[p] = 1;
	}

	// preprocess the reference sequence
	// for a bit faster access
	//UNUSED+// ssize_t r;
	ssize_t i, j ;
	for (i = 0; i < refLength; i++) {
		ref.seq[i] = nucToIndex[ref.seq[i]];
	}
	DNASequence readRC, readOrig, readRCOrig;

	ssize_t *forwardAlignment, *reverseAlignment;
	forwardAlignment = new ssize_t[read.length];
	reverseAlignment = new ssize_t[read.length];

	ssize_t numReads = 0;

	char *genomeStr;
	genomeStr= new char[readLength+1];
	genomeStr[readLength] = '\0';

  IntMatrix score;
  IntMatrix path;
  CreateMatrix(score, read.length+1, ref.length+1);
  CreateMatrix(path,  read.length+1, ref.length+1);

	do {
		if (numReads % 100 == 99)
			std::cerr << "." << std::flush;
		if (numReads % 5000 == 4999)
			std::cerr << std::endl;
		numReads++;
		unsigned char* refSeq = ref.seq;
		unsigned char* readSeq = read.seq;
		unsigned char* readRCSeq;
		ssize_t completeRead = 1;
		/*    if (read.length != readLength) {
			//			std::cout << "incomplete read\n" << std::endl;
      continue;
			}*/
		readOrig = read;
		std::string alignStr, seqStr, genomeStr;
		alignStr.resize(read.length);
		seqStr.resize(read.length);
		genomeStr.resize(read.length);
		
		MakeRC(readOrig, readRC);		
		readRCOrig = readRC;
		ssize_t nMatch, nMisMatch, nQryGap, nRefGap;
		if (gappedAlign) {
			ssize_t forScore, revScore;
			
			forScore = 
				FitAlign(read, ref, -1, 1, 5, forwardAlignment, scoreMat.scoreMat, score, path);
			revScore = 
				FitAlign(readRC, ref, -1, 1, 5, reverseAlignment, scoreMat.scoreMat, score, path);
			
			if (forScore < revScore) {
				
				std::cout << "forward: " << forScore << std::endl;
				PrintAlignment(read, ref, 0, 0, forwardAlignment, read.length, std::cout);
				ComputeAlignmentStatistics(forwardAlignment, read.length, 
																	 (char*) read.seq,
																	 (char*) ref.seq, nMatch, nMisMatch, nQryGap, nRefGap);
			}
			else {
				std::cout << "reverse: " << revScore << std::endl;
				PrintAlignment(readRC, ref, 0, 0, reverseAlignment, read.length, std::cout);
				ComputeAlignmentStatistics(reverseAlignment, read.length, 
																	 (char*) readRC.seq,
																	 (char*) ref.seq, nMatch, nMisMatch, nQryGap, nRefGap);
			}
			if (printErrorProf) {
				errProfOut << nMatch << ", " << nMisMatch << ", " << nQryGap << ", " << nRefGap << std::endl;
			}
		}
		else {

		for (j = 0; j < read.length; j++) {
			if (unmasked_nuc_index[readSeq[j]]> 3)
				completeRead = 0;
		}
		if (!completeRead)
			continue;

		MakeRC(read, readRC);
		// preprocess the read and rc for fast mismatch lookup
		readRCSeq = readRC.seq;
		for (j = 0; j < read.length; j++) {
			readSeq[j]   = nucToIndex[readSeq[j]];
			readRCSeq[j] = nucToIndex[readRCSeq[j]];
		}

		// Compute highest scoring alignment.
		double maxScore;
		ssize_t maxScorePos;
		maxScore = 0;
		maxScorePos = -1;
		double score;
		double rcScore;
		double maxRCScore = 0;
		ssize_t maxRCScorePos = -1;
		double b;
		ssize_t r;
		for (i = 0; i < refLength -read.length; i++) {
			score = 0;
			rcScore = 0;
			for (j = 0; j < read.length; j++ ) {
				// Cache some values
				r = refSeq[i+j];
				b = bias[j];
				if (r == readSeq[j])
					score+= 1;
				else 
					score -= b;
				
				if (r == readRCSeq[j])
					rcScore += 1;
				else
					rcScore -= b;
				if (-score > (read.length - j) and
						-rcScore > (read.length - j))
						break;
			}
			if (score > maxScore) {
				maxScore = score;
				maxScorePos = i;
			}
			if (rcScore > maxRCScore) {
				maxRCScore = rcScore;
				maxRCScorePos = i;
			}
			if (floor(maxScore) == read.length or 
					floor(maxRCScore) == read.length)
				break;

		}

		// print summary of alignment
		unsigned char* optSeq;
		ssize_t strand;
		ssize_t optPos;
		ssize_t numMisMatch = 0;
		strand = -1;
		if (maxScore > maxRCScore and maxScorePos >= 0) {
			optSeq = readSeq;
			strand = 0;
			optPos = maxScorePos;
			seqStr.assign((const char*) readOrig.seq, read.length);
			//			memcpy(seqStr, (const char*) readOrig.seq, read.length);
		}
		else if (maxRCScorePos >= 0) {
			optSeq = readRCSeq;
			strand = 1;
			optPos = maxRCScorePos;
			seqStr.assign((const char*) readRCOrig.seq, read.length);
			//			memcpy(seqStr, (const char*) readRCOrig.seq, read.length);
		}
		if (strand < 0) {
			continue;
		}

		std::cout << optPos << "\t" << strand << "\t";
    ssize_t index;
		for (j = 0; j < read.length; j++) {
      if (strand == 0) {
				index = j;
				genomeStr[j] = nuc_char[refSeq[optPos + j]];
			}
      else {
				index = read.length - j - 1;
				genomeStr[read.length - j - 1] = comp_ascii[nuc_char[refSeq[optPos+j]]];
			}

			if (refSeq[optPos+j] == optSeq[j]) {
        alignStr[index] = '0';
			}
			else {
  			alignStr[index] = '1';
				numMisMatch++;
			}
		}
		alignStr[read.length] = '\0';
		seqStr[read.length] = '\0';
		genomeStr[read.length] = '\0';
		std::cout << alignStr << " " << numMisMatch << "\t" 
							<< seqStr << "\t" << genomeStr << std::endl;
		}
	} 
	//while(seqReader.GetSeq(read));
  while (SeqReader::GetSeq(readsIn, read, SeqReader::noConvert));
}
		
