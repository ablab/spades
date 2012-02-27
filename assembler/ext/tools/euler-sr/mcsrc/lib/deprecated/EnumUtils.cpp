/***************************************************************************
 * Title:          EnumUtils.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "EnumUtils.h"

#include "StripGen.h"
#include "TupleLib.h"
#include "SeqUtils.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <stdlib.h>


void PrintAlignment(std::ostream *output,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    ssize_t *alignment,
		    ssize_t sign,
		    ssize_t length) {
  T_Alignment alignStruct;
  if (length == 0)
    return;
  alignStruct.locations = alignment;
  alignStruct.enumerations = &sign;
  alignStruct.length    = length;
  PrintAlignment(output, refSeq, qrySeq, alignStruct);
}

void PrintAlignment(std::ostream *output,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    T_Alignment &alignment) {

  ssize_t *enumerations, *refLocations, *qryLocations;
  if (refSeq.length == 0 || qrySeq.length == 0)
    return;
  enumerations = new ssize_t[refSeq.length];
  refLocations = new ssize_t[refSeq.length];
  qryLocations = new ssize_t[qrySeq.length];
  ssize_t numEnumerations;
  if (alignment.length < 0)
    return;
  numEnumerations =  AlignmentToEnumeration(refSeq,
					    alignment.locations,
					    Sign(alignment.enumerations[0]),
					    enumerations,
					    refLocations,
					    qryLocations);

  PrintAlignment(output, refSeq, qrySeq, enumerations, refLocations, qryLocations, numEnumerations);

  delete[] enumerations;
  delete[] refLocations;
  delete[] qryLocations;
}
		    


void PrintAlignment(std::ostream *output,
		    DNASequence &refSeq,
		    DNASequence &qrySeq,
		    ssize_t *enumerations, ssize_t *refLocations, ssize_t *qryLocations,
		    ssize_t length) {
  if (length == 0)
    return;

  ssize_t pos;
  ssize_t stripPos;
  ssize_t dir;
  pos = 0;
  // Print the alignment;
  AlignString alignmentString;
  alignmentString.lineLength = 60;
  alignmentString.outstream = output;
  
  char qryNuc, refNuc, alignChar;
  ssize_t startRefPos, startQryPos;
  
  alignmentString.startRefPos = refLocations[0];
  alignmentString.startQryPos = qryLocations[0];
  alignmentString.dir = Sign(enumerations[0]);
  ssize_t refDiff, qryDiff, minDiff;
  ssize_t refPos, qryPos;
  ssize_t refLength, qryLength;
  ssize_t i;
  pos = 0;
  dir = Sign(enumerations[pos]);
  GetAlignedLength(enumerations, refLocations, qryLocations, length, 0, refLength, qryLength);
  *output << "Enumeration : " << enumerations[0] << std::endl;
  *output << "Query length: " << qryLength 
	  << " ( " << qryLocations[0] << " ... " 
	  << qryLocations[0] + qryLength *dir + 1*dir << ")" << std::endl;
  *output << "Sbjct length: " << refLength 
	  << " ( " << refLocations[0] << " ... " 
	  << refLocations[0] + refLength << ")" <<  std::endl;

  for (pos = 0; pos < length; pos++) {
    // Start this position on two aligned nucleotides
    refPos = refLocations[pos];
    qryPos = qryLocations[pos];
    dir = Sign(enumerations[pos]);

    // Add the alignment for the position specified in the 
    // enumerations file.
    refNuc = refSeq[refPos];
    qryNuc = GetQryChar(qrySeq[qryPos], dir);
    if ((pos < length-1) && (enumerations[pos] + 1 == enumerations[pos+1])) {
      alignChar = GetAlignChar(qryNuc, refNuc);
      alignmentString.AddAlignChars(refNuc, alignChar, qryNuc, refPos, qryPos );

      // Fill in the gaps of aligned sequence until the next position in 
      // the enumeration file
      if (pos < length-1) {
	refPos++;
	qryPos+= dir;
	refDiff = refLocations[pos+1] - refPos;
	qryDiff = abs(qryLocations[pos+1] - qryPos);
	minDiff = std::min(refDiff, qryDiff);

	// Add unaligned ungapped (mismatched) characters
	for (i = 0; i < minDiff; i++, refPos++, qryPos+=dir)
	  alignmentString.AddAlignChars(tolower(refSeq.seq[refPos]), 
					' ', 
					tolower(GetQryChar(qrySeq.seq[qryPos], dir)),
					refPos, qryPos);

	// If ref pos had extra sequence, output that
	for (i = minDiff; i < refDiff; i++, refPos++)
	  alignmentString.AddAlignChars(tolower(refSeq.seq[refPos]), ' ', '-', 
					refPos, qryPos);

	// If qryPos has extra sequence, output that
	for (i = minDiff; i < qryDiff; i++, qryPos++)
	  alignmentString.AddAlignChars('-', ' ', 
					tolower(GetQryChar(qrySeq.seq[qryPos],dir)), 
					refPos, qryPos);
      }
    }
    else {
      // Done creating output for this strip.
      if (pos < length - 1) {
	// Print the strip
	alignmentString.FlushOutput(0,0);
	*output << "............................................................" 
		<<std::endl;
	// Start a new one
	alignmentString.startRefPos = refLocations[pos+1];
	alignmentString.startQryPos = qryLocations[pos+1];
	alignmentString.dir = Sign(enumerations[pos+1]);
	GetAlignedLength(enumerations, refLocations, qryLocations, 
			 length, pos+1, refLength, qryLength);
	if (refLength > 0) {
	  dir = Sign(enumerations[pos+1]);
	  *output << "Enumeration : " << enumerations[pos+1] << std::endl;
	  *output << "Query length: " << qryLength << " ( " 
		  << qryLocations[pos+1] << " ... " 
		  << qryLocations[pos+1] + qryLength *dir  << ")" << std::endl;
	  *output << "Sbjct length: " << refLength << " ( " 
		  << refLocations[pos+1] << " ... " 
		  << refLocations[pos+1] + refLength << ")" <<  std::endl;
	}
      }
      else {
	alignChar = GetAlignChar(qryNuc, refNuc);
	alignmentString.AddAlignChars(refNuc, alignChar, qryNuc, refPos, qryPos );
      }
    }
  }
  if (alignmentString.refString.length() > 0) 
    alignmentString.FlushOutput(0,0);
}


void AlignString::AddAlignChars(char refChar, char alignChar, char qryChar, ssize_t refPos, ssize_t qryPos) { 
  if (lineLength == refString.length()) 
    FlushOutput(refPos, qryPos);
  if (DNASequence::UnmaskedValue(refChar) < 4) 
    refChar = DNASequence::GetNucChar(refChar);

  if (DNASequence::UnmaskedValue(qryChar) < 4) 
    qryChar = DNASequence::GetNucChar(qryChar);
      
  refString += refChar;
  qryString += qryChar;
  alignString += alignChar;
}

void AlignString::FlushOutput(ssize_t newRefPos, ssize_t newQryPos) {
  // Print the alignment if a whole string is filled
  ssize_t i;
  *outstream << "Query:" << std::setw(7) << startQryPos << "  " << qryString;
  for (i = qryString.length(); i < lineLength; i++)
    *outstream << ' ';
  *outstream << "  " << std::setw(7) << startQryPos + dir* ((ssize_t)qryString.length()) - 1*dir << std::endl;
  *outstream << "               " << alignString;
  for (i = alignString.length(); i < lineLength; i++) *outstream << ' ';
  *outstream << "  " << std::endl;

  *outstream << "Sbjct:" << std::setw(7) << startRefPos << "  " << refString;
  for (i = refString.length(); i < lineLength; i++) *outstream << ' ';
  *outstream << "  " << std::setw(7) << startRefPos + refString.length() - 1 <<  std::endl;
  
  // Extra space between this line and the next
  *outstream << std::endl << std::endl;
  // Reset everything
  


  refString = "";
  qryString = "";
  alignString = "";
  
  
  startRefPos = newRefPos;
  startQryPos = newQryPos;
}

void GetAlignedLength(ssize_t *enumerations, ssize_t *refLocations, ssize_t *qryLocations, ssize_t length, ssize_t pos, 
		      ssize_t &refLength, ssize_t &qryLength) {
  ssize_t refStart, qryStart;
  refStart = refLocations[pos];
  qryStart = qryLocations[pos];
  while (pos < length-1) {
    if (enumerations[pos] +1 != enumerations[pos+1] || (pos+1 == length-1)) {
      refLength = abs(refLocations[pos] - refStart);
      qryLength = abs(qryLocations[pos] - qryStart);
      return;
    }
    pos++;
  }
  refLength = -1;
  qryLength = -1;
}

char GetAlignChar(char a, char b) {
  if (a >= -12 && a <= 3) {
    a = DNASequence::GetNucChar(a);
  }
  if (b >= -12 && b <= 3) {
    b = DNASequence::GetNucChar(b);
  }
  if (toupper(a) == toupper(b))
    return '|';
  else
    return ' ';
}

char GetQryChar(char nuc, ssize_t dir) {
  if (dir == 1)
    return nuc;
  else if (nuc >= -12 && nuc <= 3)
    return comp_bin[DNASequence::UnmaskedValue(nuc)];
  else
    return comp_ascii[nuc];
}

ssize_t Sign(ssize_t value) {
  if (value != 0) {
    return value / abs(value);
  }
}

ssize_t AlignmentToEnumeration(DNASequence &seq,
			   ssize_t *locations,
			   ssize_t sign,
			   ssize_t *enumerations,
			   ssize_t *refLocations,
			   ssize_t *qryLocations) {
  ssize_t i, e;
  ssize_t a; 
  a = 0;
  for (i = 0; i < seq.length; i++)  (locations[i] != -1) && a++;
  e = 0;
  for (i = 0; i < seq.length; i++) {
    if (locations[i] != -1) {
      if (sign > 0)
	enumerations[e] = (e+1);
      else
	enumerations[e] = e - a - 1;
      refLocations[e] = i;
      qryLocations[e] = locations[i];
      e++;
    }
  }
  return e;
}
