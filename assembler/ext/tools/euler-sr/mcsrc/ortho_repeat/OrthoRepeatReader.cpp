/***************************************************************************
 * Title:          OrthoRepeatReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "OrthoRepeatReader.h"

#include <iosfwd>
#include <sstream>

ssize_t OrthoRepeatReader::GetNextRepeat(OrthoRepeat *&repeat) {
  // make sure the input is ok
  if (!infile) 
    return 0;

  char refChar, qryChar, c;
  ssize_t refStrand, qryStrand;
  ssize_t refStart, refEnd, qryStart, qryEnd;
  std::string refType, qryType;

  std::stringbuf sbuf;
  // The two types of lines to parse are:
  // h chr1:1416917-1417211 AluSx c chr1:1314181-1314475 AluSp
  // h chr1:1069394-1069642 AluJb c chr1:- 


  // read in the reference half
  if (!(infile >> c)) return 0;
  if (!infile.get(sbuf, ':') or !infile.get()) return 0;
  if (!(infile >> refStart >> c >> refEnd >> refType)) return 0;


  // read in the query half
  if (!(infile >> c)) return 0;
  if (!infile.get(sbuf, ':') or !infile.get()) return 0;
  
  // The repeat may be missing from query
  if (infile.peek() == '-') {
    qryStart = qryEnd = qryStrand = -1;
    qryType = "";
  }
  else {
    if (!(infile >> qryStart >> c >> qryEnd >> qryType))
      return 0;
  }
  infile.get();
  // Parsed the entire line, create a new repeat and return
  // successful rsult
  repeat = new OrthoRepeat;
  repeat->refStart = refStart;
  repeat->refEnd   = refEnd;
  repeat->refStrand   = 0; // guess later on
  repeat->refType  = refType;
  repeat->qryStart = qryStart;
  repeat->qryEnd   = qryEnd;
  repeat->refStrand   = 0; // guess later on
  repeat->qryType  = qryType;

  return 1;
}
