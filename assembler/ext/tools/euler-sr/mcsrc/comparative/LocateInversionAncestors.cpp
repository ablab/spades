/***************************************************************************
 * Title:          LocateInversionAncestors.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <sstream>
#include <string>
#include <vector>


#include <set>

#include "CharTree.h"
#include "tree/NewettTree.h"
#include "utils.h"
#include "InversionChars.h"
#include "InversionBins.h"

int main(int argc, char* argv[]) {
  std::string inversionFile, speciesFile, treeFileName, charFileName, graphFileName;

  if (argc != 5) {
    std::cout << "usage: lia inversionFile speciesFile treeFileName graphFileName" << std::endl;
    exit(1);
  }

  inversionFile = argv[1];
  speciesFile = argv[2];
  treeFileName = argv[3];
  graphFileName = argv[4];

  std::vector<std::string> species;
  std::vector<ssize_t*> chars;
  ssize_t numChars;
  std::vector<ssize_t> invStartPositions;
  std::vector<ssize_t> invEndPositions;

  std::set<std::string> allSpecies;
  std::set<std::string> with, without;
  std::set<std::string> plus, minus;
  std::vector<std::vector<std::string>* > invSpecies;
  ssize_t i, j;
  std::vector<ssize_t> startPos, endPos;
  InversionList inversions;
  ReadInversionFile(inversionFile, species, inversions);
  ReadSpeciesFile(speciesFile, species);
  InversionList consensus;
  std::vector<ssize_t> mult, length;
  
  /*
    ReadCharsFile(charFileName, species,
    numChars, invStartPositions,invEndPositions, 
    chars, invSpecies);
  */
  CharTree cTree;
  std::ifstream treeIn;
  std::set<std::string> setA, setB;
  std::vector<std::string> *invSpecP;
  openck(treeFileName, treeIn);
  ReadTree<CharNode>(treeIn, cTree); 

  allSpecies.insert(species.begin(), species.end());
  
  // This stores a list of all leaves below each internal vertex
  cTree.StoreInternalList();
  ssize_t charNum = 0;
  ssize_t mat;
  ssize_t locus;
  for (mat = 0; mat < inversions.size(); mat++ ) {
    for (locus = 0; locus < inversions[mat]->size(); locus++) {
      // determine which species have the inversion (or differ)
      
      ssize_t found;
      plus.clear();
      minus.clear();
    
      std::set<std::string>::iterator withIt, withEnd;
      for (j = 0; j < inversions[mat]->numSpecies(); j++) {
	if (inversions[mat]->loci[locus]->orient[j] == 1) {
	  plus.insert(species[j]);
	}
	else if (inversions[mat]->loci[locus]->orient[j] == 0) {
	  minus.insert(species[j]);
	}
      }
      std::cout << "searching for charset ";
      inversions[mat]->loci[locus]->Print(std::cout);
      // Try to locate the clade in the tree for which the inversion happened.
      ssize_t minScore = 100000;
      CharNode *minCladeNode = NULL;
      CharNode *inversionNode;

      // Now locate the branch within the clade that the inversion happened on.
  
      minScore = 1000000;
      std::set<std::string>::iterator it, end;
      if (plus.size() > 0 and minus.size() > 0) {
	if (plus.size() < minus.size()) { 
	  cTree.LocateClade(plus, static_cast<CharNode*>(cTree.root), inversionNode, minScore);
	}
	else {
	  cTree.LocateClade(minus, static_cast<CharNode*>(cTree.root), inversionNode, minScore);
	}
	inversionNode->locations.push_back(inversions[mat]->loci[locus]->species);
      }
      ++charNum;
    }
  }
  std::ofstream out;
  openck(graphFileName, out, std::ios::out); 
  cTree.Print(out);
  out.close();
  return 0;
}

