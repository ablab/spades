/***************************************************************************
 * Title:          ShuffleAlign.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SHUFFLE_ALIGN_H_
#define SHUFFLE_ALIGN_H_


#include <vector>
#include "DNASequence.h"

void ShuffleAlign(DNASequence &seqA, DNASequence &seqB,
		  ssize_t nblocks,
		  ssize_t niter,
		  std::vector<double> &scores);

#endif
