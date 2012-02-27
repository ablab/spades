/***************************************************************************
 * Title:          VoteUtils.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef VOTE_UTILS_H_
#define VOTE_UTILS_H_
#include "DNASequence.h"
#include "mctypes.h"
 
ssize_t PrepareSequence(DNASequence &read);
void InitVotingMatrix(DNASequence &read, IntMatrix &votes);
void InitSolidVector(DNASequence &read, IntVector &solid);

#endif
