/***************************************************************************
 * Title:          SWAlign.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include "DNASequence.h"
#include "SeqReader.h"

#include "swat.h"



int main(int argc, char* argv[]) {
	std::cout << "this isn't finished!" << std::endl;
	assert(0);
	char *mcsrc = getenv("home");
	if (mcsrc == NULL ) {
		std::cout << "ERROR, you must set the HOME directory to where " << std::endl;
		std::cout << "projects/software/phrap/BLOSUM50 is located. " << std::endl;
	}
	

  DNASequence refSeq, qrySeq;

  if (argc != 3) {
    std::cout << "usage: swalign seq1 seq2 " << std::endl;
    return 1;
  }
  SeqReader::GetSeq(argv[0], refSeq, SeqReader::noConvert);
  SeqReader::GetSeq(argv[1], qrySeq, SeqReader::noConvert);


  Profile *q_profile;

  q_profile = make_profile_from_seq(refSeq.seq, 1, refSeq.length);
  return 0;
}
