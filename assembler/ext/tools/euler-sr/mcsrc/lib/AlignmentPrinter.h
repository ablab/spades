/***************************************************************************
 * Title:          AlignmentPrinter.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef ALIGNMENT_PRITNER_H_
#define ALIGNMENT_PRITNER_H_
#include <string>
#include "DNASequence.h"

void MakeAlignmentString(DNASequence &refSeq, 
			 DNASequence &qrySeq,
			 ssize_t refStartPos, ssize_t qryStartPos,
			 ssize_t *alignment, ssize_t alignmentLength,
			 std::string &refAlignString, 
			 std::string &qryAlignString, 
			 std::string &alignString,
			 ssize_t &refAlignStart, ssize_t &qryAlignStart);

void PrintAlignmentString(std::string refAlignString, 
			  std::string qryAlignString,
			  std::string alignString,
			  ssize_t refPos, ssize_t qryPos,
			  std::ostream &out,
			  ssize_t lineLength);


void PrintAlignment(DNASequence &refSeq, 
		    DNASequence &qrySeq,
		    ssize_t refStartPos, ssize_t qryStartPos,
		    ssize_t *alignment, ssize_t alignmentLength,
		    std::ostream &out, ssize_t lineLength=50 );

#endif
