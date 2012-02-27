/***************************************************************************
 * Title:          VoteUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "VoteUtils.h"
 
void InitVotingMatrix(DNASequence &read, IntMatrix &votes) {
	if (votes.size() < read.length) {
		CreateMatrix(votes, read.length, 9);
	}
	else {
		//UNUSED+// ssize_t j;
		ssize_t i ;
		for (i = 0; i < votes.size(); i++ ){
			std::fill(votes[i].begin(), votes[i].end(), 0);
		}
	}
}

void InitSolidVector(DNASequence &read, IntVector &solid) {
	if (read.length > solid.size()) {
		solid.resize(read.length);
	}
	std::fill(solid.begin(), solid.end(), 0);
}

ssize_t PrepareSequence(DNASequence &read) {
	ssize_t p;
	for (p = 0; p < read.length; p++ ){ 
		read.seq[p] = toupper(read.seq[p]);
		if (!(read.seq[p] == 'A' ||
					read.seq[p] == 'C' ||
					read.seq[p] == 'T' || 
					read.seq[p] == 'G'))
			return 0;
	}
	return 1;
}
