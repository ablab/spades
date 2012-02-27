/***************************************************************************
 * Title:          PrintMap.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  04/03/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "parseblast/BlastParser.h"
#include "utils.h"
#include <map>
#include <iostream>
#include <fstream>

int main(int argc, char* argv[]) {


	std::string blastName, readLengthName, mapName;
	
	if (argc != 4) {
		std::cout << "Print a map of query contigs to a reference genome." << std::endl;
		std::cout << "Usage: printMap blastFile readLengthFile mapName " << std::endl;
		exit(1);
	}
	blastName      = argv[1];
	readLengthName = argv[2];
	mapName        = argv[3];

	std::map<std::string, ssize_t> readLengths;
	std::ofstream mapOut;
	openck(mapName, mapOut, std::ios::out);

	// First read the blast file;
	BlastResult blastResult;
	ReadBlastTable(blastName, blastResult);
	
	// Now read the read lengths
	std::ifstream lenIn;
	openck(readLengthName, lenIn, std::ios::in);
	std::string readName;
	ssize_t readLength;
	while(lenIn) {
		lenIn >> readLength >> readName;
		readLengths[readName] = readLength;
	}


	// Now parse the blast file
	ssize_t b, h;
	ssize_t minMatchLength;
	ssize_t uniquelyMatched = 0;
	ssize_t multiplyMatched = 0;
	ssize_t notMatched      = 0;
	ssize_t totalHits = 0;
	ssize_t numShort = 0;
	ssize_t numDivergent = 0;
	for (b = 0; b < blastResult.size(); b++ ) {
		BlastQueryMatch *match = blastResult[b];
		if (readLengths.find(match->queryName) == readLengths.end()) {
			/*
			std::cout << "ERROR: found a match that is not in the query file " << std::endl;
			std::cout << match->queryName << std::endl;
			*/
		}
		readLength = readLengths[match->queryName];
		minMatchLength = (ssize_t) (readLength * 0.97);
		ssize_t numMatched = 0;
		for (h = 0; h < match->hsps.size(); h++ ) {
			totalHits ++;
			if (match->hsps[h].qryEnd - match->hsps[h].qryPos + 1 > minMatchLength) {
				if (match->hsps[h].identity > 0.97) {
					//UNUSED// ssize_t start;
					ssize_t refBegin, refEnd;
					if (match->hsps[h].refPos < match->hsps[h].refEnd) {
						refBegin = match->hsps[h].refPos;
						refEnd   = match->hsps[h].refEnd;
					}
					else {
					refBegin = match->hsps[h].refEnd;
					refEnd   = match->hsps[h].refPos;
					}
					
					mapOut << match->queryName << " " << refBegin - match->hsps[h].qryPos 
								 << " " << refBegin - match->hsps[h].qryPos + readLength << std::endl;
					numMatched++;
				}
				else {
					numDivergent++;
					/*
						std::cout << match->queryName << " " << match->hsps[h].identity << std::endl;
					*/
				}
			}
			else {
				/*
					std::cout << match->queryName << " " 
					<<  match->hsps[h].qryEnd - match->hsps[h].qryPos + 1 << " < " << minMatchLength << std::endl;
				*/

				numShort++;
			}
		}
		if (numMatched == 0) notMatched++;
		else if (numMatched == 1) uniquelyMatched++;
		else multiplyMatched++;
	}

	std::cout << "Processed " << totalHits << " blast hits " << std::endl;
	std::cout << uniquelyMatched << " unique hits " << std::endl;
	std::cout << multiplyMatched << " multiple hits " << std::endl;
	std::cout << notMatched << " did not have a high scoring match " << numShort << " were short "
						<< numDivergent << " were divergent " << std::endl;
	return 0;
}

