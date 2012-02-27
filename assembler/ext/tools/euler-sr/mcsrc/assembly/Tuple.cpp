/***************************************************************************
 * Title:          Tuple.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "Tuple.h"

int Tuple::tupleSize = 0;

void MutateTuple(char* tuple, ssize_t pos, unsigned char nuc) {
	tuple[pos] = nuc;
}
void DeleteTuple(char* tuple, ssize_t tupleLen, ssize_t pos, unsigned char fill) {
	ssize_t i;
	for (i = pos; i < tupleLen -1; i++ ){
		tuple[i] = tuple[i+1];
	}
	tuple[i] = fill;
}

void InsertTuple(char* tuple, ssize_t tupleLen, ssize_t pos, unsigned char nuc) {
	ssize_t i;
	for (i = tupleLen - 1; i > pos; i--) {
		tuple[i] = tuple[i-1];
	}
	tuple[i] = nuc;
}




