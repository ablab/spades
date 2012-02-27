/***************************************************************************
 * Title:          FourGameteUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef FOUR_GAMETE_UTILS_H_
#define FOUR_GAMETE_UTILS_H_

#include <map>
#include <vector>
#include <string>

#include "mctypes.h"
#include "ValidatedInversion.h"
#include "InversionBins.h"

typedef std::map<std::string, ssize_t> CharMap;
typedef std::vector<ssize_t> PosVect;
typedef std::set<ssize_t> IntSet;
typedef std::vector<IntSet > SetVector;

/*
  ssize_t FindFourGameteViolatingCharacters(IntMatrix & charMat,  
				      IntVector &inv1, IntVector &inv2);
*/

ssize_t FindFourGameteViolatingCharacters(IntMatrix & charMat,
				      IntVector & conflictCounts,
				      SetVector & conflictGraph,
				      StringVector & titles);

void GetListTitles(InversionList &inversions, StringVector &titles);
void GetListTitles(InversionMatrix &invMat, StringVector &titles);

void BuildCharMat( InversionMatrix &invMat, IntMatrix &charMat);

ssize_t RemoveMinimumConflicts(InversionMatrix &invMatrix, ssize_t maxNumConflicts);
void RemoveMinimumConflicts(InversionList &inversions, 
			    StringVector &species, 
			    ssize_t maxNumConflicts);

void GetConflictingChars(InversionMatrix &invMat, 
			 ssize_t maxNumConflicts, 
			 IntVector &conflictingChars,
			 StringVector &titles);


#endif
