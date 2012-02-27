/***************************************************************************
 * Title:          ReadPos.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/02/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/

#include "ReadPos.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "compatibility.h"

ssize_t ReadPos::hashLength;
SimpleSequenceList* ReadPos::sequences;

void PrintTuple(SimpleSequenceList &sequences, 
								ReadPos &readPosition, int tupleSize, std::ostream &out) {
  DNASequence tmpSeq;
  tmpSeq.seq = &(sequences[readPosition.read].seq[readPosition.pos]);
  tmpSeq.length = tupleSize;
  tmpSeq.PrintSeq(out);
}

/*

void PrintTuple(ssize_t tuple, int tupleSize, std::ostream &out) {
  _INT_ i;
  std::string tupleStr = "";
  _UINT_ mask = 0x3;
  unsigned char c;
  unsigned char nuc;
  for (i = 0; i < tupleSize; i++) {
    nuc = tuple & mask;
    c = nuc_char[nuc];
    tupleStr = c + tupleStr;
    tuple  >>= 2;
  }
  out << tupleStr;
}

*/
