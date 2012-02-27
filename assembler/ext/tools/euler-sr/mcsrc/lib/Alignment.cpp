/***************************************************************************
 * Title:          Alignment.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/01/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "SeqUtils.h"

void ComputeAlignmentStatistics(ssize_t *alignment, ssize_t alignLength, char *qry, char* ref,
																ssize_t &nMatch, ssize_t &nMismatch, ssize_t &nGapQry, ssize_t &nGapRef) {

	nMatch = nMismatch = nGapQry = nGapRef = 0;
	ssize_t beg, end;
	// Find the beginning of the alignment
	for (beg = 0; beg < alignLength; beg++) {
		if (alignment[beg] != -1) 
			break;
	}
	for (end = alignLength-1; end > beg; end++) {
		if (alignment[end-1] != -1)
			break;
	}

	ssize_t a;
	for (a = beg; a < end; a++ ) {
		if (alignment[a] != -1) {
			(numeric_nuc_index[(unsigned char) qry[a]] == numeric_nuc_index[(unsigned char) ref[alignment[a]]]) ? nMatch++ : nMismatch++;
			if ((a < end-1) and alignment[a] +1 != alignment[a+1]) {
				nGapRef += alignment[a+1] - alignment[a] - 1;
			}
		}
		else {
			nGapQry++;
		}
	}
}
