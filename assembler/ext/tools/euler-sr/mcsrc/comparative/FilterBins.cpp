/***************************************************************************
 * Title:          FilterBins.cpp 
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
#include <vector>
#include <set>
#include <map>
#include <sstream>

#include "utils.h"
#include "InversionBins.h"
#include "InversionUtils.h"
#include "FourGameteUtils.h"

void RemoveAllZeroRows(InversionList &invList);

int main(int argc, char*argv[]) {

  std::string seqName, inversionFileName, outFileName;
  std::string speciesFileName;
  if (argc <= 2) {
    std::cout << "usage: filterbins binIn binOut [-uo4cm]" << std::endl;
    std::cout << "-n removes null (no 0's or 1's) " << std::endl;
    std::cout << "-u removes uninformative (only 1 0 or 1 (not implemented) " 
	      << std::endl;
    std::cout << "-s species  - specify the species file name " << std::endl;
    std::cout << "-c prints the consensus sequence " << std::endl;
    std::cout << "-o orders output according to a sequence file " << std::endl;
    std::cout << "-m removes low support (deprecated) " << std::endl;
    std::cout << "-4 performs a four-gamete check on all characters " << std::endl;
    exit(1);
  }
  IntVector uninformativeCount;
  inversionFileName = argv[1];
  outFileName = argv[2];
  ssize_t printConsensus = 0;
  ssize_t removeLowSupport = 0;
  ssize_t lowSupport = 0;
  int argi = 3;
  ssize_t doFourGameteCheck = 0;
  ssize_t removeUninformative = 0;
  ssize_t printOrdered = 0;
  ssize_t keepColumn = 0;
  ssize_t removeNull = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-n") == 0) {
      removeNull = 1;
    }
    if (strcmp(argv[argi], "-u") == 0) {
      removeUninformative = 1;
    }
    if (strcmp(argv[argi], "-o") == 0) 
      printOrdered = 1;
    if (strcmp(argv[argi], "-4") == 0) 
      doFourGameteCheck = 1;
    if (strcmp(argv[argi], "-c") == 0) 
      printConsensus = 1; 
    if (strcmp(argv[argi], "-keepcol") == 0)
      keepColumn = 1;
    if (strcmp(argv[argi], "-m") == 0) {
      removeLowSupport = 1;
      lowSupport = atoi(argv[argi+1]);
      ++argi;
    }
    if (strcmp(argv[argi], "-s") == 0) {
      speciesFileName = argv[++argi];
    }
    ++argi;
  }
      
  StringVector species;
  InversionList inversions;
  ReadInversionFile(inversionFileName, species, inversions);
  if (speciesFileName != "") {
    ReadSpeciesFile(speciesFileName, species);
  }
  InversionList original;
  if (removeNull) {
    RemoveNullRows(inversions);
  }
  if (removeUninformative) {
    CopyInversionList(inversions, original);
    ssize_t i;
    for (i = 0; i < inversions.size(); i++) {
      RemoveUninformativeRows(*inversions[i], *original[i], uninformativeCount);
    }
  }
  if (doFourGameteCheck) {
    if (inversions.size() > 1) 
      RemoveMinimumConflicts(inversions, species, 0);
    else 
      RemoveMinimumConflicts(*inversions[0], 0);
  }

  std::ofstream out;
  openck(outFileName, out, std::ios::out);
  ssize_t k;
  if (!printConsensus) {
    PrintBins(species, inversions, out, printOrdered);
  }
  else {
    InversionMatrix consensus;
    GetBinConsensus(inversions, consensus);
    std::ofstream charMatFile;
    openck(outFileName, charMatFile, std::ios::out);
    PrintBins(consensus, charMatFile, species, 0);
    charMatFile.close();
  }
  out.close();
  return 0;
}

void RemoveAllZeroRows(InversionList &invList) {
  InversionList::iterator invIt;
  invIt = invList.begin(); 
  ssize_t qrySpec;
  ssize_t n1;
  ssize_t mat, inv, orient;
  for (mat = 0; mat < invList.size(); mat++ ){
    for (inv = 0; inv < invList[mat]->size(); ) {
      n1 =0 ;
      for (orient = 0; orient < invList[mat]->numSpecies(); orient++) {
	if (invList[mat]->loci[inv]->orient[orient])
	  ++n1;
      }
      if (n1 == 0) {
	invList[mat]->loci.erase(invList[mat]->loci.begin() + inv);
      }
      else {
	++inv;
      }
    }
  }
}
