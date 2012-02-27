/***************************************************************************
 * Title:          InversionBins.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "InversionBins.h"
#include <iostream>
#include "utils.h"
#include <sstream>

void ListToInvMatrix(InversionList &invList,
		     InversionMatrix &invMat) {
  ssize_t m, i;
  ValidatedInversion *inv;
  for (m = 0; m < invList.size(); m++ ) {
    for (i = 0; i < invList[m]->size(); i++ ) {
      inv = new ValidatedInversion(invList[m]->loci[i]);
      invMat.loci.push_back(inv);
    }
  }
}


void RemoveNullRows(InversionList &inversions){
  ssize_t locus, orient;
  ssize_t mat;
  for (mat = 0; mat < inversions.size(); mat++ ){
    for (locus = inversions[mat]->size() - 1; locus >= 0; --locus) {
      ssize_t num1, num0;
      num1 = num0 = 0;
      for (orient = 0; orient < inversions[mat]->numSpecies(); orient++) {
	if (inversions[mat]->loci[locus]->orient[orient] == 1) ++num1;
	if (inversions[mat]->loci[locus]->orient[orient] == 0) ++num0;
      }
      if (num1 == 0 or num0 == 0) {
	inversions[mat]->Erase(locus);
      }
    }
  }
}


void RemoveUninformativeRows(InversionMatrix &inversions,
			     InversionMatrix &original, IntVector &count) {
  std::vector<ssize_t>::iterator startIt, endIt;
  std::ofstream uninformativeOut;
  ssize_t searchChar, charNumber;
  ssize_t i;
  std::vector<ssize_t> n;
  charNumber = 0;
  //  std::cout << "RUR: starting " << std::endl;
  ssize_t locus;
  ssize_t species;
  count.resize(inversions.numSpecies());
  for (species = 0; species < count.size(); species++) {
    count[species] = 0;
  }
  for (locus = 0; locus < inversions.size(); ) {
    CountNs(inversions.loci[locus], n);
    // old bad if 
    //    if (n[1] == 1 or n[0] == 1 or (n[1] == 0 and n[0] == 0)) {
    if (n[1] <= 1 or n[0] <= 1) {
      if ((n[1] == 1 and n[0] > 0)  or (n[0] == 1 and n[1] > 0)) {
	if (n[1] == 1) {
	  searchChar = 1;
	}
	else if (n[0] == 1) {
	  searchChar = 0;
	}
	ssize_t locusIndex;
	for (locusIndex = 0; locusIndex < original.size(); locusIndex++ ) {
	  if (original.loci[locusIndex]->number == inversions.loci[locus]->number)
	    break;
	}
	assert(locusIndex < original.size());
	AnullCharacter(original.loci[locusIndex]);
	// old bad if
	if (n[1] == 1 or n[0] == 1) {
	ssize_t colIndex;
	// find wehre this was in the original list
	ssize_t specIndex;
	for (specIndex = 0; specIndex < inversions.numSpecies(); specIndex++ ) {
	  if (inversions.loci[locus]->orient[specIndex] == searchChar) {
	    break;
	  }
	}
	assert(specIndex < inversions.numSpecies());
	count[specIndex]++;
		      }
	// Only one species has this inversion, remove it.
	//      std::cout << " removing " << inversions.loci[locus]->number << std::endl;
      }
      inversions.Erase(locus);
    }
    else {
      ++locus;
    }
  }
}


void AnullCharacter(ValidatedInversion *valInv) {
  // Find the majority character
  ssize_t n0, n1;
  n0 = 0; n1 = 0;
  ssize_t i;
  for (i = 0; i < valInv->size(); i++ ) {
    if (valInv->orient[i] == 0 )
      n0++;
    if (valInv->orient[i] == 1 )
      n1++;
  }
  ssize_t majority, minority;
  if (n0 > n1) {
    majority = 0;
    minority = 1;
  }
  else {
    majority = 1;
    minority = 0;
  }
  for ( i = 0; i < valInv->size(); i++ ) { 
    if (valInv->orient[i]  == minority)
      valInv->orient[i] = majority;
  }
}


void OneMajorCharacter(ValidatedInversion *inv) {
  std::vector<ssize_t> n;
  ssize_t i;
  CountNs(inv, n);
  if (n[1] / double(n[0] + n[1]) > 0.5) {
    // Leave this alone, maybe print a horray?
  }
  else if (n[0] / double(n[0] + n[1]) > 0.5) {
    // Found something in the reverse orientation
    for (i = 0; i < inv->orient.size(); ++i) {
      if (inv->orient[i] == 0)
	inv->orient[i] = 1;
      else if (inv->orient[i] == 1)
	inv->orient[i] = 0;
      // Leave 2 alone.
    }
  }
}

void GetBinConsensus(InversionList   &inversions,
		     InversionMatrix &consensus) {
  ssize_t bin, inv;
  std::vector<ssize_t> n, n0, n1, n2;
  ssize_t invOrient;
  ssize_t i;
  ssize_t mat, locus, orient;
  if (inversions.size() == 0) {
    return;
  }
  ssize_t size = inversions[0]->numSpecies();
  for (mat = 0; mat < inversions.size(); mat++ ) {
    for (locus = 0 ; locus < inversions[mat]->loci.size(); locus++) {
      // Try to find the starting position of this inversion.
      OneMajorCharacter(inversions[mat]->loci[locus]);
    }
  }
  ssize_t nidx;
  bin = 0;
  for (mat = 0; mat < inversions.size(); mat++) {

    // Attemp to find a consensus row in this matrix
    n0.resize(size);
    n1.resize(size);
    n2.resize(size);
    for (nidx = 0; nidx < size; nidx++) {
      n0[nidx] = 0;
      n1[nidx] = 0;
      n2[nidx] = 0;
    }
    for (locus = 0; locus < inversions[mat]->size(); locus++ ) {
      for (inv = 0; inv < inversions[mat]->loci.size(); inv++) {
	invOrient = inversions[mat]->loci[locus]->orient[inv];
	if (invOrient == 0) 
	  n0[inv]++;
	else if (invOrient == 1)
	  n1[inv]++;
	else if (invOrient == 2)
	  n2[inv]++;
      }
    }
    ValidatedInversion *valInv = new ValidatedInversion;
    valInv->species = "consensus";
    valInv->orient.resize(size);

    // Try to stre the consensus if there is a clear one
    //
    ssize_t clearConsensus = 1;
    for (inv = 0; inv < size; inv++) {
      if (n0[inv] > n1[inv])
	valInv->orient[inv] = 0;
      else if (n1[inv] > n0[inv])
	valInv->orient[inv] = 1;
      else {
	clearConsensus = 0;
	valInv->orient[inv] = 2;
      }
    }

    // clear consensus implied that it is clear which
    // orientation was the consensus
    // If it wasn't clear, find the row int the inversion 
    // matrix that contains the most values, and use that
    // as the consensus
    if (clearConsensus == 0) {
      ssize_t maxFilledRow = 0;
      ssize_t maxFilled = 0;
      ssize_t numFilled;
      ssize_t row = 0;
      for (locus = 0; locus < inversions[mat]->size(); locus++) {
	numFilled = 0;
	for (inv = 0; inv < inversions[mat]->loci[locus]->size(); inv++) {
	  invOrient = inversions[mat]->loci[locus]->orient[inv];
	  if (invOrient == 0 or 
	      invOrient == 1) {
	    numFilled++;
	  }
	}
	if (numFilled > maxFilled) {
	  maxFilled = numFilled;
	  maxFilledRow = row;
	}
	++row;
      }
      for (inv = 0; inv <  inversions[mat]->loci[maxFilledRow]->size(); inv++ ){
	valInv->orient[inv] = inversions[mat]->loci[maxFilledRow]->orient[inv];
      }
    }
    consensus.loci.push_back(valInv);
  }
}

void ReadInversionFile(std::string &fileName, StringVector &species, 
		       InversionList &invList) {
  std::string token, line;
  ssize_t  start, end, bin;
  ssize_t numBins = 0;
  ssize_t size;
  std::ifstream in;
  openck(fileName, in);
  if (PeekInput(in, "species:")) {
    ReadSpeciesLine(in, species);
  }

  ValidatedInversion *inversion;
  InversionMatrix *invMatrix;
  inversion = NULL;
  invMatrix = NULL;
  ssize_t inversionNumber = 0;
  ssize_t matStarted = 0;
  while (in and in.peek() != '\n' and in.peek() != EOF) {
    if (matStarted or PeekInput(in, std::string(">"))) { 
      invMatrix = new InversionMatrix;
      invList.push_back(invMatrix);
      in >> invMatrix->title;
      std::getline(in, line);
      matStarted = 0;
    }
    else {
      if (invMatrix == NULL) {
	// Initialize the inversion if this is an anonymous inv mat
	invMatrix = new InversionMatrix;
	invList.push_back(invMatrix);
	inversionNumber = 0;
      }
      while (in) {
	if (PeekInput(in, std::string(">")) != 0) {
	  matStarted = 1;
	  break;
	}
	// read the species (each validated inverison is one locus
	// compared againsta ll other loci, so the inverison
	// should be stored.cc
	std::string species;
	if (!(in >> species)) break;
	inversion = new ValidatedInversion;
	inversion->species = species;
	// Read in the orientations of the loci in other species
	std::getline(in, line);
	if (in.eof()) return;
	std::vector<std::string> strValues;
	Tokenize(line, strValues);
	ssize_t i;
	inversion->orient.resize(strValues.size());
	inversion->number  = inversionNumber++;
	for (i = 0; i < strValues.size(); i++) {
	  if (!Convert(strValues[i], inversion->orient[i])) {
	    std::cout << "error reading inversion state " << i 
		      << "got " << inversion->orient[i] << " from " << strValues[i] 
		      << " for " << invMatrix->title << std::endl;
	  }
	  // Some backwards-compatibility
	  if (inversion->orient[i] > 2) 
	    inversion->orient[i] = 2;
	}
	invMatrix->loci.push_back(inversion);
	// sanity check
	if (invMatrix->loci.size() > 1) {
	  if (invMatrix->loci[invMatrix->loci.size()-1]->size() != 
	      invMatrix->loci[invMatrix->loci.size()-2]->size()) {
	    std::cout << "FORMAT ERROR, should have the same number of orientations per line " << std::endl;
	    std::cout << "problem with line numbers : " 
		      << invMatrix->loci[invMatrix->loci.size()-1]->size() << " "
		      << invMatrix->loci[invMatrix->loci.size()-2]->size() << std::endl;
	    std::cout << line << std::endl;
	    exit(0);
	  }
	}
	// End sanity check

      }
    }
  }
}

void PrintInversionFile(std::string &fileName, InversionList &invList) {
  std::ofstream toOut;
  openck(fileName, toOut, std::ios::out);
  ssize_t inv;
  for (inv = 0; inv < invList.size(); inv++){ 
    invList[inv]->Print(toOut);
  }
}


void PrintBins(InversionMatrix &invMatrix, 
	       std::ostream &out, StringVector &species, ssize_t ordered) {
  if (ordered) 
    invMatrix.PrintOrdered(out, species);
  else
    invMatrix.Print(out);
}

void PrintBins(StringVector &species,
	       InversionList &inversions,
	       std::ostream &out,
	       ssize_t ordered) {
  // print the species
  ssize_t s;
  for (s = 0 ; s < (ssize_t) species.size() - 1; s++) {
    out << species[s] << " ";
  }
  if (species.size() > 0)
    out << species[s];
  out << std::endl;

  // Print each bin
  ssize_t bin, inv;
  // iterating over a map, dangerous.
  ssize_t mat;
  for  (mat = 0; mat < inversions.size(); mat++) {
    PrintBins(*inversions[mat], out, species, 0);
  }
}


void CountNs(ValidatedInversion *valInv, std::vector<ssize_t> &indices, std::vector<ssize_t> &n) {
  ssize_t colIndex;
  ssize_t col;
  ssize_t invOrient;
  for (col = 0; col < indices.size(); ++col) {
    colIndex = indices[col];
    invOrient = valInv->orient[colIndex];
    if (invOrient == 0) 
      ++n[0];
    else if (invOrient == 1)
      ++n[1];
    else if (invOrient == 2)
      ++n[2];
    else if (invOrient == 3)
      ++n[3];
  }  
}


void CountNs(ValidatedInversion *valInv, std::vector<ssize_t> &n) {

  ssize_t inv, invOrient;
  n.resize(4);
  n[0] = n[1] = n[2] = n[3] = 0;
  for (inv = 0; inv < valInv->size(); ++inv) {
    invOrient = valInv->orient[inv];
    if (invOrient == 0) 
      ++n[0];
    else if (invOrient == 1)
      ++n[1];
    else if (invOrient == 2)
      ++n[2];
    else if (invOrient == 3)
      ++n[3];
  }  
  // Sanity check there should be only 0,1,2
  assert(n[0] + n[1] + n[2] + n[3] == valInv->size() or 
	 printf("%d %d %d %d %d\n", n[0], n[1], n[2], n[3], valInv->size())== 0);
}

ssize_t FindSpecies(StringVector &speciesList, std::string &species) {
  StringVector::iterator sit;
  sit = std::find(speciesList.begin(), speciesList.end(), species);
  if (sit == speciesList.end())
    return -1;
  else 
    return sit - speciesList.begin();
}

ssize_t FindSpecies(InversionMatrix &invMat, std::string &species) {
  ssize_t pos = 0;
  ssize_t i;
  for (i = 0; i < invMat.size(); i++) {
    if (invMat.loci[i]->species == species) {
      return pos;
    }
    ++pos;
  }
  return -1;
}



void IntersectLists(InversionMatrix &mata,
		    InversionMatrix &matb) {

  ssize_t a, b;
  ssize_t found;
  for (a = 0; a < mata.size(); ) {
    found = 0;
    for (b = 0; b < matb.size() and !found; b++) {
      if (mata.loci[a]->number == matb.loci[b]->number)
	found = 1;
    }
    if (!found) {
      mata.Erase(a);
    }
    else {
      a++;
    }
  }
  for (b = 0; b < matb.size(); ) {
    found = 0;
    for (a = 0; a < mata.size() and !found; a++) {
      if (matb.loci[b]->number == mata.loci[a]->number)
	found = 1;
    }
    if (!found) {
      matb.Erase(b);
    }
    else {
      b++;
    }
  }
}


void ReadSpeciesFile(std::string &specFileName, StringVector &species ) {
  std::ifstream in;
  openck(specFileName, in, std::ios::in);
  ReadSpeciesLine(in, species);
}
    

void ReadSpeciesLine(std::ifstream &in, StringVector &species) {
  std::string tempstr;
  while (in >> tempstr) {
    species.push_back(tempstr);
  }
}
