/***************************************************************************
 * Title:          SpeciesInversions.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "InversionBins.h"
#include "ValidatedInversion.h"
#include "mctypes.h"

#include <string>

void AppendSpecies(InversionList &from, InversionList &to, ssize_t fromSpecies);
void AllocateMatrices(InversionList &from, InversionList &to, ssize_t startIndex, std::string &toTitle);


int main(int argc, char* argv[]) {
  
  std::string invFileName, invOutFileName,
 specFromFileName, specToFileName;
  if (argc < 5) {
    std::cout << "usage: specinv invfile speciesfrom speciesto invout [-startindex start] [-title title]" << std::endl;
    std::cout << "       make the columns in 'invfile' coorrespond to the order " << std::endl
	<< "              and composition of the speciesto list " << std::endl;
    exit(0);
  }

  invFileName = argv[1];
  specFromFileName = argv[2];
  specToFileName   = argv[3];
  invOutFileName   = argv[4];
  ssize_t startIndex;
  std::string toTitle;
  int argi =5;
  startIndex = -1;
  toTitle    = "";
  while (argi < argc) {
    if (strcmp(argv[argi], "-startindex") == 0) {
      startIndex = atoi(argv[++argi]);
    }
    else if (strcmp(argv[argi], "-title") == 0) {
      toTitle = argv[++argi];
    }
    ++argi;
  }
  
  StringVector fromSpecies, toSpecies;
  InversionList fromInversions, toInversions;
  
  ReadInversionFile(invFileName, fromSpecies, fromInversions);
  ReadSpeciesFile(specFromFileName, fromSpecies);
  ReadSpeciesFile(specToFileName, toSpecies);
  AllocateMatrices(fromInversions, toInversions, startIndex, toTitle);
  
  ssize_t from, to;
  ssize_t specIndex;
  for (to = 0; to < toSpecies.size(); to++) {
    specIndex = -1;
    for (from = 0;from < fromSpecies.size() and specIndex == -1; from++) {
      if (fromSpecies[from] == toSpecies[to]) 
	specIndex = from;
    }
    AppendSpecies(fromInversions, toInversions, specIndex);
  }
  PrintInversionFile(invOutFileName, toInversions);
  return 0;
}

void AllocateMatrices(InversionList &from, InversionList &to, ssize_t startIndex, std::string &toTitle) {
  ssize_t inv, mat;
  for (inv = 0; inv < from.size(); inv++) {
    to.push_back(new InversionMatrix);
    for (mat = 0; mat < from[inv]->size(); mat++) {
      to[inv]->loci.push_back(new ValidatedInversion);
      if (startIndex != -1) {
	std::stringstream titleStrm;
	titleStrm << startIndex++;
	to[inv]->loci[mat]->species = titleStrm.str();
      }
      else {
	to[inv]->loci[mat]->species = from[inv]->loci[mat]->species;
      }
      to[inv]->loci[mat]->number = from[inv]->loci[mat]->number;
      if (toTitle != "") {
	to[inv]->title = toTitle;
      }
      else {
	to[inv]->title = from[inv]->title;
      }
    }
  }
}
  
void AppendSpecies(InversionList &from, InversionList &to, ssize_t species) {
  ssize_t inv, mat;
  for (inv = 0; inv < from.size(); inv++) {
    for (mat = 0; mat < from[inv]->size(); mat++) {
      if (species != -1) {
	to[inv]->loci[mat]->orient.push_back(from[inv]->loci[mat]->orient[species]);
      }
      else
	to[inv]->loci[mat]->orient.push_back(2);
    }
  }
}
