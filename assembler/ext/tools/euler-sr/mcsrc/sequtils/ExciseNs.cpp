/***************************************************************************
 * Title:          ExciseRepeats.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "utils.h"

using namespace std;

void ReplaceNs(DNASequence &seq, ssize_t minLength) {
	//UNUSED// ssize_t i = 0;
	//UNUSED// ssize_t startN = -1;
	unsigned char*end = &seq.seq[seq.length];
	unsigned char *cp, *np;
	np = NULL;
	for (cp = &seq.seq[0]; cp != end; ++cp) {
		if (*cp == 'N') {
			if (np == NULL) {
			np = cp;
			}
		}
		else {
			if (np != NULL) {
				// Found a stretch of N's.
				if (cp - np < minLength) {
					for (; np != cp ; np++) {
						*np = RandomNuc();
					}
				}
				// Not in a stretch of N's any more.
				np = NULL;
			}
		}
	}
}


int main(int argc, char *argv[]) {

  if (argc != 3) {
		cout << "usage: exciseNs infile outfile " << endl;
		cout << "   [-replace L (10)]  Don't excise stretches of 'N' that " << endl
				 << "                 are less than L long.  Instead, replace " << endl
				 << "                 them by random nucleotides." << endl;
    exit(1);
  }
  std::string inFileName, outFileName;

  inFileName = argv[1];
  outFileName= argv[2];
	int argi = 3;
	ssize_t minLength = 10;
	while (argi < argc) {
		if (strcmp(argv[argi], "-replace") == 0) {
			minLength = atoi(argv[++argi]);
		}
		++argi;
	}
  std::ifstream in;
  std::ofstream out;
  
  openck(inFileName, in);
  openck(outFileName, out, std::ios::out);

  DNASequence seq, norep;

  while (SeqReader::GetSeq(in, seq, SeqReader::noConvert)) {
    if (seq.length == 0) 
      break;

		ReplaceNs(seq, minLength);

    //UNUSED// ssize_t nonMasked = 0;
    unsigned char *cur;
    unsigned char *notN;
    unsigned char *end = &seq.seq[0] + seq.length;
		//UNUSED// ssize_t index =  0;
		//UNUSED// ssize_t length = 0;
		notN = cur = seq.seq;
		while (cur != end) {
			// Skip a masked sequence.
			while (cur != end and  
						 (*cur == 'X' or
							*cur == 'N'))
				cur++;

			// Go to end of good sequence
			while (cur != end and *cur != 'X' and *cur != 'N') {
				*notN  = *cur;
				cur++;
				notN++;
			}
		}
		seq.length = notN - seq.seq;
		seq.PrintlnSeq(out);
  }
}
  





  
  
    
