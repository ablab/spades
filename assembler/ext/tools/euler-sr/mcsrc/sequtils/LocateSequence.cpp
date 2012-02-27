/***************************************************************************
 * Title:          LocateSequence.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/05/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SimpleSequence.h"
#include "SeqReader.h"
#include <string>
using namespace std;
int main(int argc, char* argv[]) {
	DNASequenceList seqList;
	if (argc < 3) {
		cout << "usage: locateSeq string sub" << endl;
		exit(1);
	}

	std::string seqFileName, seq;
	seqFileName = argv[1];
	seq = argv[2];

	ReadDNASequences(seqFileName, seqList);

	ssize_t s;
	ssize_t p;
	ssize_t c;
	ssize_t compLen = seq.size();
	const char *seqPtr = seq.c_str();
	ssize_t match = 0;
	for (s = 0; s < seqList.size(); s++ ) {
		for (p = 0; p < seqList[s].length - seq.size() + 1; p++ ) {
			match = 1;
			for (c = 0; c < compLen; c++ ) {
				if (seqList[s].seq[p+c] != seqPtr[c]) {
					match = 0;
					break;
				}
			}
			if (match) {
				std::cout << seqList[s].namestr << " " << s << " " << p << std::endl;
			}
		}
	}
	return 0;
}
				
