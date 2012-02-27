/***************************************************************************
 * Title:          CountWords.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/10/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include <vector>
#include <string>

using namespace std;
ssize_t IsMasked(unsigned char* c, ssize_t len) {
	unsigned char *end;
	end = c + len;

	while (c != end) {
		if (*c == 'N' or *c == 'X') return 1;
		c++;
	}
	return 0;
}

class SortWords {
public:
	ssize_t wordLength;
	unsigned char *str;
	ssize_t operator()(ssize_t a, ssize_t b) {
		return strncmp((const char*) &str[a], (const char*) &str[b], wordLength) < 0;
	}
};
int main(int argc, char* argv[]) {
	string inName= argv[1];
	ssize_t wordLen = atoi(argv[2]);

	vector<ssize_t> words;

	DNASequence chr;
	SeqReader::GetSeq(inName, chr, SeqReader::noConvert);
	unsigned char *c, *end;
	end = &chr.seq[chr.length];
	c = chr.seq;
	while (c != end) { *c = toupper(*c); ++c;}
	end = &chr.seq[chr.length - wordLen ];
	ssize_t numWords = chr.length - wordLen ;
	words.resize(numWords);
	ssize_t i;
	i = 0;
	for (i = 0; i < numWords;i++) {
		words[i] = i;
	}
		
	SortWords sorter;
	sorter.wordLength = wordLen;
	sorter.str = chr.seq;
	cout << "sorting." << endl;
	std::sort(words.begin(), words.end(), sorter);

	// Count the # of unique words.
	ssize_t i1, i2;;
	i1 =0 ;i2 = 0;
	ssize_t numUnique = 0;
	while (i1 < numWords and i2 < numWords) {
		if (IsMasked(&chr.seq[words[i1]], wordLen)) {
			i1++;
			continue;
		}
		i2 = i1;
		while (strncmp((const char*) &chr.seq[words[i1]], 
									 (const char*) &chr.seq[words[i2]], wordLen) == 0)
			i2++;
		numUnique++;
		i1 = i2;
	}
	cout << "counted: " << numUnique << " unique words of length: " << wordLen << endl;
}
	
