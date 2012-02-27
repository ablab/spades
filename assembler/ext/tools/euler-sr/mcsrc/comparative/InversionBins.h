/***************************************************************************
 * Title:          InversionBins.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _INVERSION_BINS_H_
#define _INVERSION_BINS_H_

#include "ValidatedInversion.h"
#include <set>

ssize_t FindSpecies(InversionMatrix &invMat, std::string &species);

ssize_t FindSpecies(StringVector &speciesList, std::string &species);

void CountNs(ValidatedInversion *valInv, std::vector<ssize_t> &n); 

void CountNs(ValidatedInversion *valiInv, std::vector<ssize_t> &indices, std::vector<ssize_t> &n);


void ListToInvMatrix(InversionList &invList,
		     InversionMatrix &invMat);

void GetBinConsensus(InversionList   &inversions,
		     InversionMatrix &consensus);

void ReadInversionFile(std::string &fileName, StringVector &species, 
		       InversionList &invList);

void PrintInversionFile(std::string &fileName, InversionList &invList);

void PrintBins(InversionMatrix &invMatrix, 
	       std::ostream &out, StringVector &species, ssize_t ordered = 0);

void PrintBins(StringVector &species,
	       InversionList &inversions,
	       std::ostream &out,
	       ssize_t ordered);

void ReadSpeciesLine(std::ifstream &in, StringVector &species);
void ReadSpeciesFile(std::string &specFileName, StringVector &species );

void RemoveNullRows(InversionList &inversions);

void RemoveUninformativeRows(InversionMatrix &inversions,
			     InversionMatrix &original,
			     IntVector &species);

void IntersectLists(InversionMatrix &mata,
		    InversionMatrix &matb);

void AnullCharacter(ValidatedInversion *valInv);


#endif
