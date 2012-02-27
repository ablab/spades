/***************************************************************************
 * Title:          ThreadUtils.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  04/13/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef THREAD_UTILS_H_
#define THREAD_UTILS_H_

#include "IntervalGraph.h"
#include "ThreadPath.h"
#include "mctypes.h"
#include <string>
#include <vector>

ssize_t CollectBranchSequences(IntervalGraph &graph, ssize_t branchLength, ssize_t vertex, ssize_t maxBranches,
													 std::vector<char*> &branchSequences, 
													 std::vector<ssize_t> &branchLengths, 
													 string curSeq, ssize_t &curSeqIndex, ssize_t depth);


ssize_t StorePathSequences(IntervalGraph &graph, 
											 ssize_t edgeIndex,ssize_t edgePos, 
											 std::string curSequence,ssize_t searchLength, 
											 std::vector<std::string> &sequences);

ssize_t Thread(IntervalGraph &graph,ssize_t edgeIndex, ssize_t intvIndex, 
					 const char* read, ssize_t readLength,
					 ThreadPath &minThreadPath,
					 ssize_t maxScore, ssize_t maxThreadLength, 
					 DNASequence &perfectRead,
					 std::vector<ssize_t> &scoreList, ssize_t scoreListLength, IntMatrix &scorMat,
					 ssize_t readIndex, ssize_t allowGaps);

void AppendReverseComplements(DNASequenceList &sequences);

void FindBestMatch(DNASequence &read, std::vector<std::string> &pathSequences, 
									 ssize_t & id_best, ssize_t & delta_best, ssize_t &delta_best2);

ssize_t ThreadRecursively(IntervalGraph &graph,
											ssize_t edgeIndex, ssize_t edgePos,
											const char*curSeq,ssize_t curSeqLength,
											ThreadPath &curThreadPath, ssize_t curThreadScore,
											ThreadPath &minThreadPath, ssize_t &minThreadScore, 
											ssize_t maxScore, ssize_t maxDepth, DNASequence &perfectRead, 
											std::vector<ssize_t> &scoreList, ssize_t scoreListLength, ssize_t depth, 
											IntMatrix &scoreMat, ssize_t printNext, ssize_t readIndex, ssize_t allowGaps);


ssize_t ThreadRecursivelyBandedAlign(DNASequence &read,
																 IntervalGraph &graph,
																 ssize_t edgeIndex, ssize_t edgePos,
																 char* curPathSeq, ssize_t curPathSeqLength,
																 ThreadPath &curThreadPath, ssize_t curThreadScore,
																 ThreadPath &minThreadPath, ssize_t minThreadScore,
																 ssize_t curDepth, ssize_t maxDepth, map<ssize_t,ssize_t> &nTraverse, string spacing);


void ThreadToSeq(IntervalGraph &graph,
								 ThreadPath &path,
								 std::string &sequence);

#endif
