/***************************************************************************
 * Title:          ExtractLoci.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include "DNASequence.h"
#include "FragmentUtils.h"
#include "utils.h"

int main(int argc, char* argv[] ) {
  
  std::string locusFileName;
  std::string sequenceDir;
  std::string outFileName;
  if (argc < 3) { 
    std::cout << "usage: extractloci locusFile locusSeqfile [-seqdir"
              << " sequencedir] " << std::endl;
    return 1;
  }  
  int argi= 1;
  locusFileName = argv[argi++];
  outFileName   = argv[argi++];

  sequenceDir = ".";
  while (argi < argc) {
    if (strcmp(argv[argi], "-seqdir") == 0) {
      ++argi;
      sequenceDir = argv[argi];
    }
    ++argi;
  }

  Loci loci;

  ReadLoci(locusFileName, loci);
  ReadLociSequences(loci, sequenceDir);

  std::ofstream locusOut;
  openck(outFileName, locusOut, std::ios::out);
  ssize_t i;
  for (i = 0; i < loci.size(); i++) {
    loci[i]->locusSeq.PrintSeq(locusOut);
    locusOut << std::endl;
  }
  return 0;
}
