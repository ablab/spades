/***************************************************************************
 * Title:          EnumUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _ENUM_UTILS_
#define _ENUM_UTILS_

#include "DNASequence.h"
#include "Alignment.h"

#include <iostream>
#include <fstream>
#include <string>

class AlignString {
public:
  std::ostream *outstream;
  std::string refString;
  std::string qryString;
  std::string alignString;
  ssize_t startRefPos;
  ssize_t startQryPos;
  ssize_t lineLength;
  ssize_t dir;
  void FlushOutput(ssize_t newRefPos, ssize_t newQryPos);
  void AddAlignChars(char refChar, char alignChar, char qryChar, ssize_t refPos, ssize_t qryPos);
  AlignString() {refString = qryString = alignString = "";}
};

char GetQryChar(char nuc, ssize_t dir);
char GetAlignChar(char a, char b);

void GetAlignedLength(ssize_t *enumerations, ssize_t *refLocatiosn, ssize_t *qryLocations, 
		      ssize_t length, ssize_t pos, 
		      ssize_t &refLength, ssize_t &qryLength);

void PrintAlignment(std::ostream *output,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    ssize_t *alignment,
		    ssize_t sign,
		    ssize_t length);

void PrintAlignment(std::ostream *outfile,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    ssize_t *enumerations, ssize_t *refLocations, ssize_t *qryLocations,
		    ssize_t length);

void PrintAlignment(std::ostream *output,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    T_Alignment &alignment);

ssize_t AlignmentToEnumeration(DNASequence &seq,
			   ssize_t *locations,
			   ssize_t sign,
			   ssize_t *enumerations,
			   ssize_t *refLocations,
			   ssize_t *qryLocations);

ssize_t Sign(ssize_t value);
	       

#endif
