/***************************************************************************
 * Title:          ReadMap.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "ReadMap.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include <iterator>
#include <iostream>
#include "ParseTitle.h"

void ReadReadMapFile(std::string &mapFileName, 
										 ReadMapList &readMap) {

	std::ifstream in;
	openck(mapFileName, in, std::ios::in);
	ReadMap rm;
	std::string rem, posStr;
	//UNUSED// ssize_t pos;
	ssize_t cur = 0;

	while(in.good()) {
		if (in.get() != '>')
			return;
		readMap.push_back(rm);
		in >> readMap[cur].name;
		std::getline(in, rem);
		in >> readMap[cur].start >> readMap[cur].end;
		std::getline(in, rem);
		/*		std::getline(in, rem);
		std::getline(in, posStr);
		std::stringstream posStrm(posStr);
		
		copy(std::istream_iterator<ssize_t>(posStrm), std::istream_iterator<ssize_t>(),
				 back_inserter(readMap[cur].map));
		*/
		cur++;
	}
}

ssize_t ReadMappedReads(std::string readFileName, 
										ReadMapList &readMap,
										DNASequenceList &reads) {

	DNASequence read;
	std::ifstream readIn;
	std::string readName;
	openck(readFileName, readIn, std::ios::in);
	ssize_t curMapPos = 0;
	while(SeqReader::GetSeq(readIn, read, SeqReader::noConvert) and curMapPos < readMap.size()) {
		// the read map should be in the same order as the list of reads
		ParseTitle(read.namestr, readName);
		if (readName == readMap[curMapPos].name) {
			reads.push_back(read);
			curMapPos++;
		}
	}
	return reads.size();
}




