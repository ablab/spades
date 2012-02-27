/***************************************************************************
 * Title:          PairReads.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>

#include <map>
#include <string>
#include <list>
#include "utils.h"

ssize_t IsTitle(std::string &title) {
	if (title.size() == 0) return 0;
	else return (title.c_str()[0] == '>');
}

ssize_t ParsePairedTitle(std::string &title, std::string &base, std::string &dir) {
	base = "";
	dir  = "";
	if (title.size() == 0) return 0;
	base = title;
	if (base.c_str()[0] == '>') 
		base.replace(0,1,"");

	ssize_t dotPos = base.find(".");
	if (dotPos != base.npos) {
		dir  = base.substr(dotPos+1);
		base = base.substr(0, dotPos);
		return 1;
	}
	else {
		return 0;
	}
}

typedef std::map<std::string, ssize_t> posMap;
int main(int argc, char* argv[]) {

	std::string readsFile;
	std::string pairOutFile;
	ssize_t rule;
	
	int argi = 1;
	if (argc < 3) {
		std::cout << "usage: pairReads readsFile pairOutFile rule" << std::endl;
		exit(1);
	}
	readsFile = argv[argi++];
	pairOutFile = argv[argi++];
	rule = atoi(argv[argi++]);
	
	std::ifstream readsIn;
	openck(readsFile, readsIn, std::ios::in);

	std::ofstream pairsOut;
	openck(pairOutFile, pairsOut, std::ios::out);

	std::string line, base, dir;
	ssize_t readIndex = 0;
	std::list<std::string> mateTitles;
	std::list<ssize_t> mateIndices;
	posMap names;
	while (readsIn) {
		std::getline(readsIn, line);
		if (IsTitle(line)) {
			if (ParsePairedTitle(line, base, dir)) {
				if (dir == "x") {
					names[base] = readIndex;
				}
				else if (dir == "y") {
					mateTitles.push_back(base);
					mateIndices.push_back(readIndex);
				}
			}
			++readIndex;
		}
	}
	
	//UNUSED// ssize_t mate;
	std::list<std::string>::iterator titleIter;
	std::list<ssize_t>::iterator indexIter;

	for (titleIter = mateTitles.begin(),
				 indexIter = mateIndices.begin(); titleIter != mateTitles.end(); 
			 ++titleIter, ++indexIter) {
		if (names.find(*titleIter) != names.end()) {
			pairsOut << *titleIter << " " << names[*titleIter] << " " << *indexIter << " " << rule << std::endl;
		}
	}
	return 0;
}
