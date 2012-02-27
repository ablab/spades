/***************************************************************************
 * Title:          ExtractInversions.cpp 
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
#include "InversionFile.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"

void PrintUsage() {
  std::cout << "extractinversions  invFile refSpecies qrySpecies "
	    << "sequenceName invSequenceFileName [-q]" << std::endl;
}

int main(int argc, char* argv[]) {

  if (argc < 6) {
    PrintUsage();
    exit(1);
  }

  std::string invFileName;
  std::string refSpecies, qrySpecies;
  std::string sequenceName, invSequenceFileName;
  ssize_t useRef = 1;
  invFileName = argv[1];
  refSpecies  = argv[2];
  qrySpecies  = argv[3];
  sequenceName = argv[4];
  invSequenceFileName = argv[5];
  int argi = 6;
  ssize_t doHardMask = 1;
  while (argi < argc) {
    if (strcmp(argv[argi], "-q") == 0) {
      useRef = 0;
      std::cout << "using query sequence " << std::endl;
    }
    if (strcmp(argv[argi], "-n") == 0) {
      doHardMask = 0;
    }
    ++argi;
  }
  InversionMap inversions;

  StringSet  species;

  ReadInversionFile(invFileName, inversions, species);
  StringSet::iterator sit, send;
  for (sit = species.begin(); sit != species.end(); sit++ ) {
    std::cout << "got species: " << (*sit) << std::endl;
  }
  std::cout << "got " << inversions.size() << " inversions " << std::endl;
  DNASequence seq, invSeq;
  DNASequence rc;
  SeqReader::MaskRepeats();
  if (!SeqReader::GetSeq(sequenceName, seq, SeqReader::noConvert)) {
    exit(1);
  } 
  if (doHardMask) 
    seq.HardMask();
  if ( seq.length == 0) {
    std::cout << "No sequence read " << std::endl;
    exit(1);
  }
  
  if (!useRef) {
    MakeRC(seq, rc);
    seq = rc;
    rc.Free();
  }
  invSeq._ascii = 1;
  std::ofstream invSequenceFile;
  openck(invSequenceFileName, invSequenceFile, std::ios::out);
  if (inversions.find(refSpecies) != inversions.end()) {
    if (inversions[refSpecies].find(qrySpecies) != 
	inversions[refSpecies].end()) {
      InversionVector::iterator specInvIt, specInvEnd;
      specInvEnd = inversions[refSpecies][qrySpecies].end();
      specInvIt  = inversions[refSpecies][qrySpecies].begin();

      for (; specInvIt != specInvEnd; ++specInvIt) {
	if (useRef) {
	  invSeq.seq = &seq.seq[(*specInvIt)->tStart];
	  invSeq.length = (*specInvIt)->tEnd - (*specInvIt)->tStart + 1;
	  *invSeq.titlestream << refSpecies << "_" 
			     << (*specInvIt)->tStart << "_" 
			     << (*specInvIt)->tEnd << "_" 
			     << qrySpecies;
;
	}
	else {
	  invSeq.seq = &seq.seq[(*specInvIt)->qStart];
	  invSeq.length = (*specInvIt)->qEnd - (*specInvIt)->qStart + 1;
	  *invSeq.titlestream << refSpecies << "_" << qrySpecies << "_" 
			     << (*specInvIt)->qStart << "_" 
			     << (*specInvIt)->qEnd;
	}
	invSeq.PrintSeq(invSequenceFile);
	invSequenceFile << std::endl;
      }
      invSequenceFile.close();
      return 0;
   }
    std::cout << "query species "  << qrySpecies  
	      << " not found " << std::endl;
    return 0;
  }
  std::cout << "reference " << refSpecies 
	    << " not found " << std::endl;
  return 1;
}
