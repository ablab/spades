/***************************************************************************
 * Title:          SOLiD.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SOLiD.h"
#include "utils.h"

void SOLiD::ReadCSFile(std::string fileName, 
											 std::vector<std::string> &sequences,
											 std::vector<std::string> &titles) {

	std::ifstream inFile;
	openck(fileName, inFile, std::ios::in);
	std::string line;
	// Discard all lines starting with a '#'
	while(inFile.peek() == '#') {
		std::getline(inFile, line);
	}
	std::string seq, title;
	while(inFile) {
		if (!std::getline(inFile, line))
			break;
		if (line.c_str()[0] == '>') {
			title = line.substr(1);
		}
		else {
			SOLiD::CSToNuc((const char*) &(line.c_str()[1]), line.size() - 1,
										 SOLiD::NucIndex[line.c_str()[0]], seq);
			sequences.push_back(seq);
			titles.push_back(title);
		}
	}
}

ssize_t SOLiD::CSToNuc(const char* colorSeq, 
									 ssize_t colorSeqLength,
									 unsigned char firstNuc,
									 std::string &seq) {

	char *nucSeq;
	seq = "";
	nucSeq = new char[colorSeqLength+1];
	nucSeq[colorSeqLength] = '\0';
	ssize_t i;
	ssize_t curColor = firstNuc;
	ssize_t nextColor;
	char curNuc;
	curNuc = firstNuc;
	for (i = 0; i < colorSeqLength; i++ ) {
		nextColor = SOLiD::LetterToNumber[colorSeq[i]];
	  nucSeq[i] = SOLiD::SolidTransNuc[curNuc][nextColor];
		curNuc = SOLiD::NucIndex[nucSeq[i]];
	}
	seq = nucSeq;
	delete[] nucSeq;
}
