//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * read_reconstruction.cpp
 *
 *  Created on: Nov 28, 2011
 *      Author: snikolenko
 */

#include "standard.hpp"
#include "read_reconstruction.hpp"
#include "globals.hpp"
#include "hammer_tools.hpp"

using namespace std;

void ReadReconstructor::reconstruct(const string & fname) {
	TIMEDLN(" ");

}

void ReadReconstructor::reconstruct(const string & fn_left, const string & fn_right) {

}

void ReadReconstructor::reconstructAll() {
	changed_reads = 0;
	changed_bases = 0;
	int iter_no = 0;
	while (makeOneReconstructionIteration()) {
		++iter_no;
		TIMEDLN("Iteration " << iter_no << " is over. Changed " << changed_bases << " bases in " << changed_reads << " so far");
	}
	TIMEDLN("Iterative reconstruction is over");
}

bool ReadReconstructor::makeOneReconstructionIteration() {
	for (size_t i = 0; i < Globals::pr->size(); ++i) {
		if ( read_corrected[i] ) continue;
		PositionRead & r = (*Globals::pr)[i];
		read_corrected[i] = correctOneRead(r);
	}
	return (changed_reads > 10);
}

bool ReadReconstructor::correctOneRead(PositionRead & pr) {
   assert(false);
   return false;
}
