/***************************************************************************
 * Title:          SplitPairedFastq.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <fstream>
#include <stdlib.h>
#include "utils.h"

using namespace std;
int main(int argc, char* argv[]){ 
	string joinedFastqFile = argv[1];
	ssize_t    forReadLength = atoi(argv[2]);
	string splitFastqFile  = argv[3];

	ifstream fastqIn;
	ofstream fastqOut;
	openck(joinedFastqFile, fastqIn, std::ios::in);
	openck(splitFastqFile, fastqOut, std::ios::out);

	string seqTitle, qualTitle;
	string seq, qual;
	string forSeq, forQual, revSeq, revQual;
	while(fastqIn) {
		//		char c = fastqIn.get(); // get the @
		fastqIn.get(); // get the @
		// TODO: should verify it's "@"

		getline(fastqIn, seqTitle); fastqIn.get();

		getline(fastqIn, seq); fastqIn.get();

		// get the "+\n". TODO: should verify it's "+\n"
		fastqIn.get(); fastqIn.get();
		getline(fastqIn, qual); fastqIn.get();
		if (seq.size() < forReadLength*2) {
			cout << "ERROR, sequence: " << endl
					 << seq << endl
					 << " is too short." << endl;
			exit(0);
		}
		forSeq.assign(seq, 0, forReadLength);
		revSeq.assign(seq, forReadLength, forReadLength);
		forQual.assign(qual, 0, forReadLength);
		revQual.assign(qual, forReadLength, forReadLength);
		fastqOut << "@" << seqTitle << "/1" << endl;
		fastqOut << forSeq << endl;
		fastqOut << "+" << seqTitle << "/1" << endl;
		fastqOut << forQual << endl;
		fastqOut << "@" << seqTitle << "/2" << endl;
		fastqOut << revSeq << endl;
		fastqOut << "+" << seqTitle << "/2" << endl;
		fastqOut << revQual << endl;
	}
}

	
