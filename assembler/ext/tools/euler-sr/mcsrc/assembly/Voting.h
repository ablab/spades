/***************************************************************************
 * Title:          Voting.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  04/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef VOTING_H_
#define VOTING_H_

#include "SimpleSequence.h"
#include "SeqReader.h"
#include "SeqUtils.h"
#include "utils.h"
#include "hash/HashUtils.h"
#include "BufferedSeqReader.h"
#include "FixErrorsStats.h"
#include "mctypes.h"
#include "VoteUtils.h"

#include <vector>
#include <iostream>
#include <iomanip>
#include <ext/hash_map>
#include <map>
#include <deque>
#include "IntegralTuple.h"
#include "VoteUtils.h"

ssize_t FindSolidSubsequence(DNASequence &seq, CountedIntegralTupleDict &spectrum, 
												 int tupleSize, ssize_t minMult, ssize_t &seqStart, ssize_t &seqEnd);

ssize_t TrimSequence(DNASequence &seq, CountedIntegralTupleDict &spectrum,	int tupleSize, 
								 ssize_t minMult, ssize_t &seqStart, ssize_t &seqEnd, ssize_t maxTrim, ssize_t printMap);

ssize_t FixSequence(DNASequence &seq, 
								//								T_Spectrum &spectrum,
								CountedIntegralTupleDict &spectrum,
								IntMatrix &votes, IntVector &alreadySolid,
								int tupleSize, ssize_t minMult, ssize_t voteThreshold, Stats &stats, ssize_t numSearch,
								ssize_t &fixPos);

ssize_t VoteSequence(DNASequence &seq, CountedIntegralTupleDict &spectrum, int tupleSize,
								 ssize_t minMult, ssize_t startPos,
								 IntMatrix &votes, IntVector &solid, 
								 ssize_t numSearch, ssize_t runFast,
								 ssize_t checkInsertions, ssize_t checkDeletions, ssize_t earlyEnd,
								 std::deque<ssize_t> &history);

void VoteHistory(IntMatrix &votes, std::deque<ssize_t> &history);

ssize_t CheckSolid(DNASequence &seq, CountedIntegralTupleDict &spectrum, int tupleSize, ssize_t minMult);

ssize_t SolidifySequence(DNASequence &read, 
										 CountedIntegralTupleDict &spectrum, int tupleSize, ssize_t minMult,
										 IntMatrix &votes, IntVector &solid,
										 ssize_t minVotes, Stats &stats, ssize_t numSearch, 
										 ssize_t DoDeletion, ssize_t DoInsertion, ssize_t earlyEnd);

#endif
