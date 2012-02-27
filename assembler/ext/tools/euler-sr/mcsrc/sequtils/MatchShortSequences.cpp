/***************************************************************************
 * Title:          MatchShortSequences.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/10/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include <vector>
#include <string>

class Compare {
public:
	static ssize_t len;
	ssize_t operator()(char *a, char*b) {
		return strncmp(a, b, len);
	}
};

ssize_t Compare::len = 35;
int main(int argc, char* argv[]) {
	std::string genomeName, seqListName;
	genomeName = argv[1];
	seqListName = argv[2];

	DNASequence genome, genomeRC;
	SeqReader::GetSeq(genomeName, genome, SeqReader::noConvert);
	Compare::len = 35;
	std::vector<char*> ptrList;
	ssize_t i;
	ssize_t readLen = 35;
	ptrList.resize(genome.length* 2 - (readLen-1)*2);
	for (i = 0; i < ptrList.size(); i++) 
		ptrList[i] = NULL;

	MakeRC(genome, genomeRC);
	
	for (i = 0; i < genome.length - (readLen) + 1; i++) {
		ptrList[i] = (char*) &(genome.seq[i]);
	}
	std::cout << "last pos: " << i  << " genome len: " <<  genome.length << std::endl;
	for (i = 0; i < genomeRC.length - (readLen) + 1 ; i++) {
		ptrList[i + genomeRC.length - (readLen)+  1] = (char*) &(genomeRC.seq[i]);
	}
	Compare comp;
	std::cout << "sorting..." << std::endl;
	for (i = 0; i < ptrList.size(); i++) {
		assert(ptrList[i] != NULL);
		ptrList[i] = ptrList[i];
	}

	std::sort(ptrList.begin(), ptrList.end(), comp);

	std::ifstream in;
	openck(seqListName, in, std::ios::in);
	DNASequence query;
	std::vector<char*>::iterator it;
	while (SeqReader::GetSeq(in, query, SeqReader::noConvert)) {
		it = std::lower_bound(ptrList.begin(), ptrList.end(), (char*) query.seq);
		if (it != ptrList.end() and strncmp(*it, (char*)query.seq, 35)) {
			query.PrintlnSeq(std::cout);
		}
	}
}
	
	
