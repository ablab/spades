/***************************************************************************
 * Title:          AxtReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <sstream>
#include "AxtReader.h"
#include <stdlib.h>

void AxtReader::ReadAxtFile(std::string &axtFileName, 
			    AxtEntries &axtEntries) {
  std::ifstream axtIn;
  axtIn.open(axtFileName.c_str());
  if (!axtIn.good()) {
    std::cout << "could not open " << axtFileName << std::endl;
    exit(0);
  }
  AxtEntry *axtEntry;
  while (GetAxtEntry(axtIn, axtEntry)) {
    axtEntries.push_back(axtEntry);
  }
}

ssize_t AxtReader::GetAxtEntry(std::ifstream &axtIn,
			   AxtEntry *&axtEntry) {
  if (axtIn.eof()) {
    axtEntry = NULL;
    return 0;
  }
  // parse the title
  axtEntry = new AxtEntry;
  if (!AxtReader::GetTitle(axtIn, *axtEntry)) { delete axtEntry; axtEntry = NULL; return 0; }
  
  std::vector<ssize_t> refGaps, qryGaps;
  if (!GetGapLocations(axtIn, axtEntry->refGapLocations)) { delete axtEntry; axtEntry = NULL; return 0; }
  if (!GetGapLocations(axtIn, axtEntry->qryGapLocations)) { delete axtEntry; axtEntry = NULL; return 0; }

  axtIn.get(); // next line is a '\n';
  return !axtIn.eof();
}

ssize_t AxtReader::IsNuc(char c) {
  return (c == 'a' || c == 't' || c == 'g' || c == 'c' ||
	  c == 'A' || c == 'T' || c == 'G' || c == 'C');
}
	  
ssize_t AxtReader::GetGapLocations(std::ifstream &axtIn, 
			       std::vector<ssize_t> &gapLocations) {
  //UNUSED// ssize_t gaplocation = 0;
  ssize_t position    = 0;
  ssize_t prevWasDNA  = 0;
  char c = '\0';
  // Initialize the dna status
  // Want to record the transitions from gapped states to non-gapped 
  // (actg), so know the current state (c), and the state at the last 
  // transition.  
  c = axtIn.get();
  if (c == '-') 
    prevWasDNA = 0;
  else {
    prevWasDNA = 1;
    gapLocations.push_back(position);
  }
  while (c != '\n' && c != '\0' && !axtIn.eof()) {
    if (c == '-' && prevWasDNA) {
      gapLocations.push_back(position-1);
      prevWasDNA = 0;
    }
    else if (IsNuc(c) && !prevWasDNA) {
      prevWasDNA = 1;
      gapLocations.push_back(position);
    }
    c = axtIn.get();
    position++;
  }
  // If the sequence didn't end on a gap, end a block here.
  if (prevWasDNA) {
    gapLocations.push_back(position-1);
  }
  return !axtIn.eof();
}

ssize_t AxtReader::GetTitle(std::ifstream &axtIn, AxtEntry &axtEntry) {
  char strand;
  axtIn >> axtEntry.number 
	>> axtEntry.refTitle >> axtEntry.refStart >> axtEntry.refEnd
	>> axtEntry.qryTitle >> axtEntry.qryStart >> axtEntry.qryEnd
	>> strand >> axtEntry.score;
  if (strand == '+') 
    axtEntry.qryStrand = 0;
  else
    axtEntry.qryStrand = 1;

  if (axtIn.eof()) return 0;
  
  //  axtIn.get(sb, '\n');
  //  std::cout << "extracted to: " << sb.str() << std::endl;
  axtIn.get();
  return !axtIn.eof();
}
  
