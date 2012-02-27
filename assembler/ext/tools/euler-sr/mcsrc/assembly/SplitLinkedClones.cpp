/***************************************************************************
 * Title:          SplitLinkedClones.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "compatibility.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "ParseTitle.h"
#include "SeqUtils.h"
#include "align/alignutils.h"
#include "AlignmentPrinter.h"
#include "IntegralTupleStatic.h"


void PrintUsage() {
	std::cout << "usage: splitLinkedClones seqFile linkerFile pairFile" << std::endl;
	std::cout << "  -forName 'fName' ('a') Name clone forward reads with fastaName_fName" << std::endl;
	std::cout << "  -revName 'rName' ('b') Name clone reverse reads with fastaName_rName"
						<< std::endl;
	std::cout << "  -startRead 'index' (0) Start indexing read mates at 'index'" 
						<< std::endl;
	std::cout << "  -minIdentity 'I' (0.70) Only split clones if the identity is at last I%"
						<< endl << endl;
	std::cout << "  -minReadLength 'l' (40) Only extract paired reads if the length " << endl << endl
						<< "                         of each side is greater than 'l'" << endl << endl;
	std::cout << "  -singletons 'file' (none) Print singletons to 'file'" << endl << endl;
	std::cout << "  -notPal (false)        Assumes the linker sequence is not palindromic. " << endl
						<< "                         Normally, linkers such as the 454 linker are palindromes" << endl 
						<< "                         so that only one orientation of the sequene should be aligned" << endl
						<< "                         to the linker.  If this is not the case, the orientation" << endl
						<< "                         of the sequence is assumed to be the one with the highest " << endl
						<< "                         alignment score." << endl;
}


ssize_t FindLinker( DNASequence &linker, DNASequence &read, 
								ssize_t palindromicLinker, DNASequence &linkerRC,
								Score &matchScores, 
								IntMatrix &scoreMat, IntMatrix &pathMat,
								IntMatrix &rcScoreMat, IntMatrix &rcPathMat,
								ssize_t *alignment, ssize_t *revAlignment, ssize_t *&optAlignment,
								double minIdentity, ssize_t &linkerStart, ssize_t &linkerEnd) {
		
	if (scoreMat.size() < read.length) {
		CreateMatrix(scoreMat, linker.length + 1, read.length + 1);
		CreateMatrix(pathMat, linker.length + 1, read.length + 1);
		if (!palindromicLinker) {
			CreateMatrix(rcScoreMat, linker.length + 1, read.length + 1);
			CreateMatrix(rcPathMat, linker.length + 1, read.length + 1);
		}
	}
	double alignScore;
	alignScore    = FitAlign(linker, read, -1, 2, 3, alignment, matchScores.scoreMat,
													 scoreMat, pathMat);
			
	//	cout << alignScore << endl;
	if (palindromicLinker) {
		// make it so that the reverse alignment is never chosen
		optAlignment = alignment;
	}
	else {
		// Actually have to compute the alignment score to determine
		// the direction of the linker.
		double revAlignScore;
		revAlignScore = FitAlign(linkerRC, read, -1, 2, 3, revAlignment, matchScores.scoreMat,
														 rcScoreMat, rcPathMat);
				
		if (revAlignScore < alignScore) 
			optAlignment = revAlignment;
	}

	ssize_t alignBegin, alignEnd;
			
	alignBegin = 0;
	alignEnd = linker.length;
	// FInd the boundaries of the fit.
	while(alignBegin < linker.length and optAlignment[alignBegin] == -1)
		alignBegin++;
	while (alignEnd > 0 and alignEnd > alignBegin and optAlignment[alignEnd-1] == -1)
		alignEnd--;
			
			
	linkerStart = alignment[alignBegin];
	linkerEnd   = alignment[alignEnd-1];
	ssize_t linkerAlignmentLength = linkerEnd - linkerStart + 1;
	/*	cout << alignBegin << " " << alignEnd << endl;
	PrintAlignment(linker, read, 0, 0, alignment, linker.length, cout);
	cout << ((double) -alignScore)/linkerAlignmentLength << endl;
	*/
	// Look to see if there is a really long alignment

	if ( (((double) -alignScore)/linkerAlignmentLength) < minIdentity or
			 alignBegin >= alignEnd or
			 -alignScore < 10) {
		return 0;
	}
	else {
		return 1;
	}
}
int main(int argc, char* argv[]) {

	std::string readFileName, linkerFileName, pairFileName;

	std::string forCloneName, revCloneName;

	forCloneName = "a";
	revCloneName = "b";

	DNASequence read;
	DNASequence linker, linkerRC;

	if (argc < 4) {
		PrintUsage();
		exit(1);
	}
	int argi = 1;
	readFileName = argv[argi++];
	linkerFileName = argv[argi++];
	pairFileName = argv[argi++];
	double minIdentity = 0.70;
	ssize_t readNumber = 0;
	ssize_t minReadLength = 40;
	ssize_t palindromicLinker = 1;
	std::string singletonFileName = "";
	std::ofstream singletonFile;
	while (argi < argc ) {
		if (strcmp(argv[argi], "-forName") == 0){ 
			forCloneName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-revName") == 0){ 
			revCloneName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-startRead") == 0){ 
			readNumber = atosz(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-singletons") == 0) {
			singletonFileName = argv[++argi];
		}
		else if (strcmp(argv[argi], "-minReadLength") == 0) {
			minReadLength = atoi(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-minIdentity") == 0) {
			minIdentity = atof(argv[++argi]);
		}
		else if (strcmp(argv[argi], "-notPal") == 0) {
			palindromicLinker = 0;
		}
		else {
			PrintUsage();
			std::cout << "bad option: " << argv[argi] << std::endl;
			exit(1);
		}
		++argi;
	}

	if (singletonFileName != "") {
		openck(singletonFileName, singletonFile, std::ios::out);
	}

	SeqReader::GetSeq(linkerFileName, linker, SeqReader::noConvert);
	MakeRC(linker, linkerRC);

	std::ifstream readFile;
	std::ofstream readMateFile, mateFile;
	
	openck(readFileName, readFile,     std::ios::in);
	openck(pairFileName, readMateFile, std::ios::out);
	
	Score score(-1,3,1,1);
	ssize_t *alignment = NULL;
	ssize_t *revAlignment = NULL;
	ssize_t alignmentSize = 0;
	

	IntMatrix scoreMat, rcScoreMat;
	IntMatrix   pathMat, rcPathMat;
	//UNUSED// ssize_t maxReadLength = 0;
	//UNUSED// ssize_t alignScore, revAlignScore;
	DNASequence forward, reverse, reverseRC;
	std::string readName;
	
	//		if (alignmentSize < read.length) {
	//			if (alignment != NULL) {
	delete[] alignment;
	delete[] revAlignment;
				//			}
	alignment = new ssize_t[linker.length+1];
	revAlignment = new ssize_t[linker.length+1];
	alignmentSize = linker.length;
	//		}
	ssize_t readIndex = 0;
	while(SeqReader::GetSeq(readFile, read, SeqReader::noConvert)) {
		ssize_t i;
		PrintStatus(readIndex);
		readIndex++;
		for (i = 0; i < alignmentSize; i++) {
			alignment[i] = -1;
			revAlignment[i] = -1;
		}
		
		ssize_t foundLinker;
		ssize_t linkerStart, linkerEnd;
		ssize_t *optAlignment;
		//		cout << "first find linker:" << endl;
		foundLinker = FindLinker(linker, read, palindromicLinker, linkerRC,
														 score, 
														 scoreMat, pathMat, rcScoreMat, rcPathMat, 
														 alignment, revAlignment, optAlignment,
														 minIdentity, linkerStart, linkerEnd);

		readName = "";
		ParseTitle(read.namestr, readName);


		if (!foundLinker) {
			// the linker wasn't found in the sequence, just output 
			// this as a singleton.
			if (singletonFileName != "") {
				read.PrintlnSeq(singletonFile);
			}
			else {
				read.PrintlnSeq(readMateFile);
			}
		}
		
		else {
			/*
				cout << read.namestr << endl;

			*/
			// Find the beginning of the linker in the read.
			
			// If the linker is at the very end of the read, it is not
			// paired, discard it here for simplicity
			/*
			if (alignEnd <= alignBegin or
					optAlignment[alignEnd-1] == read.length-1) {
				std::cout << "the linker should not be in the end of a read but it seems to be." << std::endl;
				continue;
			}
			*/
			

			forward.seq = &read.seq[linkerEnd+1];
			forward.length = read.length - linkerEnd-1;
			
			ssize_t forLinkerFound = 0;
			do {
				ssize_t forLinkerStart, forLinkerEnd;
				//				cout << "forward search: " << endl;
				forLinkerFound = FindLinker(linker, forward, palindromicLinker, linkerRC, score,
																		scoreMat, pathMat, rcScoreMat, rcPathMat,
																		alignment, revAlignment, optAlignment,
																		minIdentity, 
																		forLinkerStart, forLinkerEnd);
				if (forLinkerFound) {
					//					cout << read.namestr << " found forward linker: " << forLinkerStart << " " << forLinkerEnd << endl;
					// find the longer side of the sequence outside the linker,
					// and use that as the read.
					if (forLinkerStart > forward.length - forLinkerEnd) {
						forward.length = forLinkerStart;
					}
					else {
						forward.seq = &forward.seq[forLinkerEnd+1];
						forward.length = forward.length - forLinkerEnd - 1;
					}
				}
			} while (forLinkerFound);

			if (forward.length > minReadLength) {
				forward.namestr = readName + "." + forCloneName;
				forward.PrintlnSeq(readMateFile);
				readNumber++;
			}


			// The first part of the read is actually the end of the clone.
			reverse.seq    = &read.seq[0];
			reverse.length = linkerStart;

			ssize_t revLinkerFound = 0;
			do {
				ssize_t revLinkerStart, revLinkerEnd;
				//				cout << "reverse search " << endl;
				revLinkerFound = FindLinker(linker, reverse, palindromicLinker, linkerRC,
																		score,
																		scoreMat, pathMat, rcScoreMat, rcPathMat,
																		alignment, revAlignment, optAlignment,
																		minIdentity, revLinkerStart, revLinkerEnd);
				if (revLinkerFound) {
					//					cout << read.namestr << " found a reverse linker " << revLinkerStart << " " << revLinkerEnd << endl;
					if (revLinkerStart > (reverse.length - revLinkerEnd - 1)) {
						reverse.length = revLinkerStart;
					}
					else {
						reverse.seq = &reverse.seq[revLinkerEnd+1];
						reverse.length = reverse.length - revLinkerEnd - 1;
					}
				}
			}
			while (revLinkerFound);

			if (reverse.length > minReadLength) {
				MakeRC(reverse, reverseRC);
				reverseRC.namestr = readName + "." + revCloneName;
				reverseRC.PrintlnSeq(readMateFile);
				readNumber++;
				reverseRC.Reset(0);
			}
		}
	}
}
