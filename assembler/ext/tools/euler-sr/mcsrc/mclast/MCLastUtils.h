/***************************************************************************
 * Title:          MCLastUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef  MCLAST_UTILS_H_
#define  MCLAST_UTILS_H_

#include "DNASequence.h"
#include "hash/HashUtils.h"
#include "MCHsp.h"
#include "mctypes.h"

ssize_t ParseHashOptions(int argc, char**argv,
		     HashPattern &hashPattern, 
		     ssize_t &mask);

void InitHashPattern(std::string fileName, HashPattern &hashPattern);

void ReadHashPatternString(std::string &hashPatternStr,
			   HashPattern &hashPattern);

void ReadHashPattern(std::string &patFileName,
		     HashPattern &hashPattern);

ssize_t FindCommonK(DNASequence &refSeq, DNASequence &qrySeq,
		HashTable &refHashTable, HashTable &qryHashTable,
		HashPattern &hashPattern, 
		std::vector<ssize_t> &refPos, std::vector<ssize_t> &qryPos,
		ssize_t maskRepeats=0);

ssize_t FindHsps(DNASequence &refSeq, DNASequence &qrySeq,
	     HashTable &refHashTable, HashTable &qryHashTable,
	     ssize_t window,
	     HashPattern &hashPattern,
	     FloatMatrix &scoreMat, double gapOpen, double gapExtend,
	     ssize_t maxScore, ssize_t maxDistance,
	     std::vector<MCHsp> &hsps);

double ExtendMatch(DNASequence &refSeq, ssize_t refPos,
		  DNASequence &qrySeq, ssize_t qryPos,
		  HashPattern &hashPattern,
		  HashValue &hashValue,
		  FloatMatrix &scoreMat, double gapOpen, double gapExtend,
		  MCHsp &hsp,
		  ssize_t *&alignment, ssize_t &alignentLength);

void InitHashPattern(HashPattern &hashPattern);

#endif
