/***************************************************************************
 * Title:          WindowGC.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqReader.h"
#include "DNASequence.h"

using namespace std;
int main(int argc, char* argv[]) {
	string seqName;
	ssize_t window;
	if (argc < 3) {
		cout << "usage: windowGC seqName window" << endl;
		exit(0);
	}
	seqName = argv[1];
	window  = atoi(argv[2]);
	DNASequence sequence;
	SeqReader::GetSeq(argv[1], sequence, SeqReader::noConvert);

	ssize_t i;
	ssize_t gc;
	char *sptr;
	for (i = 0; i < sequence.length - window + 1; i++ ){
		char* end = (char*) &sequence.seq[i] + window;
		char c;
		gc = 0;
		for (sptr = (char*) &sequence.seq[i]; sptr != end; sptr++) {
			c = *sptr;
			if (c == 'g' or 
					c == 'G' or
					c == 'c' or
					c == 'C') {
				++gc;
			}
		}
		cout.precision(2);
		cout <<  double(gc)/window << endl;
	}
}
			
