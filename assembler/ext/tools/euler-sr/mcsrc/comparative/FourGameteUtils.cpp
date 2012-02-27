/***************************************************************************
 * Title:          FourGameteUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "FourGameteUtils.h"

#include <map>
#include <set>
#include <vector>
#include "utils.h"

void GetListTitles(InversionMatrix &invMat, StringVector &titles) {
  ssize_t locus;
  for (locus = 0; locus < invMat.size(); locus++ ) {
    titles.push_back(invMat.loci[locus]->species);
  }
}

void GetListTitles(InversionList &inversions, StringVector &titles) {
  ssize_t mat;
  ssize_t locus;
  for (mat = 0; mat < inversions.size(); mat++ ) {
    GetListTitles(*inversions[mat], titles);
  }
}
  
void RemoveMinimumConflicts(InversionList &inversions, StringVector &species, ssize_t maxNumConflicts) {

  // Transform the binMap to an inversion list
  InversionList invList;
  InversionMatrix invMat;
  StringVector titles;
  ListToInvMatrix(inversions, invMat);
  GetListTitles(inversions, titles);
  //  std::cout << "searching inv list of size " << invMat.size() << std::endl;
  IntVector conflictingChars;
  
  GetConflictingChars(invMat, maxNumConflicts, conflictingChars, titles);
  //  std::cout << "got " << conflictingChars.size() << " cc's " << std::endl;

  // Remomve the conflicting chars from the bin map
  ssize_t charIndex = 0;
  ssize_t removeIndex = 0;
  // Make sure the conflicting chars are in sorted order
  // since the 
  if (conflictingChars.size() == 0) return;

  std::sort(conflictingChars.begin(), conflictingChars.end());
    
  ssize_t invIndex = 0;
  ssize_t binInv; 
  ssize_t conflictIndex = 0;
  ssize_t bin = 0;
  ssize_t charsErased;
  IntVector conflictingIndices;
  ssize_t mat, inv;
  for (mat = 0; 
       mat < inversions.size() and 
	 conflictIndex < conflictingChars.size() ; ) {
    charsErased = 0;
    conflictingIndices.clear();
    for (inv = 0; inv < inversions[mat]->size(); ) {
      if (invIndex == conflictingChars[conflictIndex]) {
	ssize_t specIndex;
	specIndex = FindSpecies(species, inversions[mat]->loci[inv]->species);
	conflictingIndices.push_back(specIndex);
	inversions[mat]->loci.erase(inversions[mat]->loci.begin() + inv);
	++conflictIndex;
	++charsErased;
      }
      else {
	++inv;
      }
      ++invIndex;
    }
    if ( inversions[mat]->size() > 0 and conflictingIndices.size() > 0 ) {
      // Some characters have been removed from the bin. (inversion locus)
      // That implies that these characters represet unreliable
      // alignments, so they should be removed from the rest of the
      // alignments in thisbin.
      ssize_t locus;
      for (locus = 0; locus < inversions[mat]->size(); locus++) {
	ssize_t conflict;
	for (conflict = 0; conflict < conflictingIndices.size(); conflict++ ) 
	  inversions[mat]->loci[locus]->orient[conflictingIndices[conflict]] = 2;
      }
    }
    
    if (inversions[mat]->size() == 0) {
      //      std::cout << "removing inversion " << mat << std::endl;
      inversions.erase(inversions.begin() + mat);
    }
    else {
      ++mat;
    }
  }
}

void GetConflictingChars(InversionMatrix &invMat, ssize_t maxNumConflicts, IntVector &conflictingChars, StringVector &titles) {
  IntMatrix charMat;
  IntVector ind1, ind2;
  IntVector counts;
  ssize_t numConflicts = 1;
  IntVector maxIndices;
  
  SetVector conflictGraph;
  IntVector conflictCount;
  IntVector toRemove;
  IntVector indices;
  conflictGraph.resize(invMat.size());
  conflictCount.resize(invMat.size());
  BuildCharMat(invMat, charMat);

  ssize_t i;
  FindFourGameteViolatingCharacters(charMat, conflictCount, conflictGraph, titles);
  // Count the number of conflicts;

  for (i = 0; i < conflictCount.size(); i++) {
    std::cout << "char: " << titles[i] << " " << conflictCount[i] << std::endl;
    if (conflictCount[i] > 0) 
      numConflicts++;
  }
  
  while (numConflicts >= maxNumConflicts) {
    ssize_t maxConflictingIndex = 0;
    ssize_t maxConflicts;
    // Find the most conflicting character.
    maxConflicts = 0;
    maxConflictingIndex = -1;
    for (i = 0; i < conflictCount.size(); i++){
      if (maxConflicts < conflictCount[i]) {
	maxConflictingIndex = i;
	maxConflicts = conflictCount[i];
      }
    }
    //    std::cout << "max cnflicts: " << maxConflicts << std::endl;
    if (maxConflictingIndex < 0)
      break;
    //    std::cout << "conflict at: " << maxConflictingIndex << std::endl;
    conflictingChars.push_back(maxConflictingIndex);
    IntSet::iterator cftIt, cftEnd;

    // Now erase every edge that goes to this vertex
    // and decrement the number of conflicts this vertex has.
    cftEnd = conflictGraph[maxConflictingIndex].end();
    cftIt  = conflictGraph[maxConflictingIndex].begin();
    std::cout << "vertex  "<< titles[maxConflictingIndex] << " " << maxConflicts << std::endl;
    for (; cftIt != cftEnd; ++cftIt) {
      // sanity checks.
      // We're not removing ourselves, and 
      // we're not removing something that doesn't exist.
      assert(*cftIt != maxConflictingIndex);
      assert(conflictGraph[*cftIt].find(maxConflictingIndex) !=
	     conflictGraph[*cftIt].end());
      
      if (maxConflicts == 1) {
	// Only two characters conflict with eachother, remove both of these
	conflictingChars.push_back(*cftIt);
	std::cout << " ... also removing "<< titles[*cftIt] << std::endl;
      }

      conflictCount[*cftIt]--;
      conflictGraph[*cftIt].erase(maxConflictingIndex);
    }
    conflictCount[maxConflictingIndex] = 0;
    numConflicts = 0;
    for (i = 0; i < conflictCount.size(); i++) {
      if (conflictCount[i] > 0) 
	numConflicts++;
    }
  }
}

ssize_t RemoveMinimumConflicts(InversionMatrix &invMatrix, ssize_t maxNumConflicts) {
  IntVector toRemove;
  StringVector titles;
  GetListTitles(invMatrix, titles);
  GetConflictingChars(invMatrix, maxNumConflicts, toRemove, titles);
  // Remove the chars from the inv list
  std::sort(toRemove.begin(), toRemove.end());
  ssize_t i;
  for (i = toRemove.size()-1; i >= 0; i--) {
    std::cout << "removing conflict " << invMatrix.loci[toRemove[i]]->species 
	      << " " << toRemove[i] << std::endl;
    invMatrix.loci[toRemove[i]]->Print();
    invMatrix.Erase(toRemove[i]);
  }
  return toRemove.size();
}


ssize_t FindFourGameteViolatingCharacters(IntMatrix & charMat,
				      IntVector & conflictCounts,
				      SetVector & conflictGraph,
				      StringVector &titles) {
  ssize_t s;
  ssize_t numInversions = charMat.size();
  if (numInversions <= 0) 
    return 0;

  ssize_t numSpecies = charMat[0].size();

  // Find species/traits that violate the 4-gamete test.

  ssize_t i1, i2;  // Inversion features

  ssize_t aa, ab, ba, bb;
  // Look through every pair of inversions.
  std::ofstream conflictGraphStream;
  openck("conflicts.dot", conflictGraphStream, std::ios::out);
  conflictGraphStream << "digraph G { " << std::endl;
  
  for (i1 = 0; i1 < numInversions-1; i1++) {
    for (i2 = i1 + 1; i2 < numInversions; i2++) {
      aa = 0; ab = 0; ba = 0; bb = 0;
      for (s = 0; s < numSpecies; s++) {
	if (charMat[i1][s] == 0 and charMat[i2][s] == 0)
	  aa = 1;         		          
	if (charMat[i1][s] == 0 and charMat[i2][s] == 1) 
	  ab = 1;         		          
	if (charMat[i1][s] == 1 and charMat[i2][s] == 0)
	  ba = 1;         		          
	if (charMat[i1][s] == 1 and charMat[i2][s] == 1)
	  bb = 1;
      }
      if (aa and ab and ba and bb) {
	// Found a violation of 4-gamete condition. How should I output this? 
	conflictCounts[i1]++;
	conflictCounts[i2]++;
	ssize_t c;
	conflictGraph[i1].insert(i2);
	conflictGraph[i2].insert(i1);

	conflictGraphStream << titles[i1] << " -> " << titles[i2] << ";" << std::endl;
      }
    }
  }
  conflictGraphStream << "}" << std::endl;
  conflictGraphStream.close();
  return 1;
}

void BuildCharMat( InversionMatrix &invMat, IntMatrix &charMat) {
  ssize_t bin, row;
  row = 0;
  ssize_t numLoci  = 0;
  ssize_t numSpecies = 0;

  numLoci  = invMat.size();
  numSpecies = invMat.numSpecies();
  CreateMatrix(charMat, numLoci, numSpecies);
  ssize_t mat, locus, spec;
  for (locus = 0; locus < numLoci; locus++ ) {
    for (spec = 0; spec < numSpecies; spec++) {
      charMat[row][spec] = invMat.loci[locus]->orient[spec];
    }
    ++row;
  }
}
