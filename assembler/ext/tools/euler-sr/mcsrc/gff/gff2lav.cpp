/***************************************************************************
 * Title:          gff2lav.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <string.h>
#include "SeqReader.h"
#include "DNASequence.h"

#include "GFFFile.h"
#include "lav/LAVFile.h"
#include "lav/LAVPrinter.h"

int main(int argc, char* argv[]) {

  std::string gffInName, lavOutName;

  GFFFile gffFile;
  
  if (argc < 3) {
    std::cout << "usage: gff2lav gffFile lavFile -rs ... -qs ... [-rl ...] [-ql ...]" << std::endl;
		std::cout << "  -rs reference_sequence_name" << std::endl;
		std::cout << "  -qs query_sequence_name" << std::endl;
		std::cout << "  -rl reference_length" << std::endl;
		std::cout << "  -ql query_length" << std::endl;
		std::cout << "       -rl and -ql are computed automatically from -rs and -qs" << std::endl;
    exit(0);
  }

  gffInName = argv[1];
  lavOutName = argv[2];

  ssize_t refLength = 0, qryLength = 0;
  int argi = 3;
  DNASequence refSeq;
  DNASequence qrySeq;
  while (argi < argc) {
    if (strcmp(argv[argi], "-rl") == 0) {
      refLength = atoi(argv[argi+1]);
      argi+=2;
    }
    else if (strcmp(argv[argi], "-ql") == 0) {
      qryLength = atoi(argv[argi+1]);
      argi+=2;
    }
    else if (strcmp(argv[argi], "-rs") == 0) {
      std::string refSeqName = argv[argi+1];
      SeqReader::GetSeq(refSeqName, refSeq);
      refLength = refSeq.length;
      argi+=2;
    }
    else if (strcmp(argv[argi], "-qs") == 0) {
      std::string qrySeqName = argv[argi+1];
      SeqReader::GetSeq(qrySeqName, qrySeq);
      qryLength = qrySeq.length;
      argi+=2;
    }
  }
  gffFile.ParseGFFFile(gffInName);

  LAVFile lavFile;
  gffFile.ToLAVFile(lavFile, 
		    refSeq.namestr, refLength, 
		    qrySeq.namestr, qryLength);
  
  std::ofstream lavOut;
  openck(lavOutName, lavOut, std::ios::out);
  LAVPrinter::PrintLAVFile(lavFile, lavOut);
  return 0;
}
