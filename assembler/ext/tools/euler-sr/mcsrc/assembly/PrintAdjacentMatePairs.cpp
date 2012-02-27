/***************************************************************************
 * Title:          PrintAdjacentMatePairs.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "MateLibrary.h"
#include "utils.h"
#include "IntegralTupleStatic.h"


int main(int argc, char* argv[]) {
	std::string readsName, mateTableName, readsOutName;
	if (argc < 3) {
		std::cout << "usage: printAdjacentMatePairs readsFile mateTable readsOut" << std::endl;
		exit(1);
	}
	
	readsName = argv[1];
	mateTableName = argv[2];
	readsOutName = argv[3];

	std::ofstream readsOut;
	openck(readsOutName, readsOut, std::ios::out);
	ReadMateList mateTable;
	ReadMateTable(mateTableName, mateTable);
	std::cout << "done reading mate table." << std::endl;
	//UNUSED// ssize_t readIndex = 0;

	//DNASequenceList reads;
	//	ReadDNASequences(readsName, reads);
	std::vector<std::string> titles;
	std::vector<std::string> reads;
	std::ifstream readsIn;
	openck(readsName, readsIn, std::ios::in);
	while (readsIn) {
		std::string title, read;
		std::getline(readsIn, title);
		std::getline(readsIn, read);
		titles.push_back(title);
		reads.push_back(read);
	}
		
	std::cout << "done reading reads." << std::endl;
	ssize_t r;
	ssize_t mateIndex;
	for (r = 0; r < reads.size(); r++ ) {
		if (reads[r].size() > 0) {
			mateIndex = mateTable[r].mateIndex;
			if (mateIndex != -1) {
				//reads[r].PrintlnSeq(readsOut);
				readsOut << titles[r] << std::endl;
				readsOut << reads[r] << std::endl;
				readsOut << titles[mateIndex] << std::endl;
				readsOut << reads[mateIndex] << std::endl;

				// Make sure this isn't double printed.
				mateTable[r].mateIndex = -1;
				mateTable[mateIndex].mateIndex = -1;
			}
		}
	}
}
