/***************************************************************************
 * Title:          RepeatUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef REPEAT_UTILS_H_
#define REPEAT_UTILS_H_

#include "DNASequence.h"
#include "ProfileCount.h"

ssize_t DetermineOrientation(DNASequence &ref, DNASequence &qry, ssize_t *&locations, double **scoreMat);
void  StoreProfile(DNASequence &consensus, DNASequence &repeatSeq, 
		   ProfileCount &profileCount, double **scoreMatPtr);
#endif
