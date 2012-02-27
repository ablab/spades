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
#include "utils.h"

int main(int argc, char *argv[]) {

  if (argc != 3) {
    std::cout << "usage: excrep infile outfile " << std::endl;
    exit(1);
  }
  std::string inFileName, outFileName;

  inFileName = argv[1];
  outFileName= argv[2];
	int argi = 3;
	ssize_t minLength = 0;
	while (argi < argc) {
		if (strcmp(argv[argi], "-minLength") == 0) {
			minLength = atoi(argv[++argi]);
		}
		++argi;
	}
  std::ifstream in;
  std::ofstream out;
  
  openck(inFileName, in);
  openck(outFileName, out, std::ios::out);

  DNASequence seq, norep;

  SeqReader::MaskRepeats();

  while (SeqReader::GetSeq(in, seq, SeqReader::noConvert)) {
    if (seq.length == 0) 
      break;
    //UNUSED// ssize_t nonMasked = 0;
    unsigned char *c;
    unsigned char *p;
    unsigned char *end = &seq.seq[0] + seq.length;
		ssize_t index =  0;
		//UNUSED// ssize_t length;
		c = seq.seq;
		while (c != end) {
			// Skip a repeat.
			while (c != end and  
						 ((*c >= 'a' and *c <= 't') or
							*c == 'X' or
							*c == 'N'))
				c++;
			p = c;
			
			// Go to end of good sequence
			while (c != end and *c >= 'A' and *c <= 'T' and *c != 'X' and *c != 'N') 
				c++;
				
			ssize_t len = c - p;
			
			norep.seq = p;
			norep.length = len;

			*norep.titlestream << seq.namestr << "_norep_" << index;
			norep.PrintSeq(out);
			out << std::endl;
			index++;
		}
  }
}
  





  
  
    
