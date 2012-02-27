/***************************************************************************
 * Title:          SRSPairParser.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <sstream>
#include "utils.h"

#include "SRSPairParser.h"

const ssize_t SRSPairParser::LENGTH = 1;
const ssize_t SRSPairParser::IDENTITY = 2;
const ssize_t SRSPairParser::SIMILARITY = 3;
const ssize_t SRSPairParser::GAPS = 4;
const ssize_t SRSPairParser::SCORE = 5;
const ssize_t SRSPairParser::COMMENT = 6;
const ssize_t SRSPairParser::BLANK = 7;
const ssize_t SRSPairParser::SEQUENCE = 8;
const ssize_t SRSPairParser::ALIGN = 9;


ssize_t SRSPairParser::GetCommentLineType(std::string &line) {
  if (line.find("Length") != line.npos) {
    return LENGTH;
  }
  else if (line.find("Identity") != line.npos) {
    return IDENTITY;
  }
  else if (line.find("Similarity") != line.npos ) {
    return SIMILARITY;
  }
  else if (line.find("Gaps") != line.npos ) {
    return GAPS;
  }
  else if (line.find("Score") != line.npos )
    return SCORE;
  else
    return 0;
}

ssize_t SRSPairParser::ParseFractionLine(std::string &line, ssize_t &num, ssize_t &denom, double &frac) {
  ssize_t cPos;
  cPos = line.find(":");
  if (cPos <= 0) 
    return 0;
  std::istringstream  lineStream(line);
  std::string temp;
  char c;
  if (! (lineStream >> temp))  return 0;
  if (! (lineStream >> temp)) return 0;
  if (! (lineStream >> num) ) return 0;
  if (! (lineStream >> c) ) return 0;
  if (! (lineStream >> denom ) ) return 0;
  if (! (lineStream.get() ) ) return 0;
  if (! (lineStream.get() ) ) return 0;
  if (! (lineStream >> frac) ) return 0;
  return 1;
}

ssize_t SRSPairParser::GetLine(std::ifstream &in, std::string &line, ssize_t repeat) {
  char buf[1000];
  ssize_t i;
  for (i = 0; i < repeat; i++) {
    if (in.getline(buf,1000) ) {
      line = buf;
    }
    else {
      line = "";
      return 0;
    }
  }
  return 1;
}

ssize_t  SRSPairParser::GetLineType(std::ifstream &in) {
  char c;
  c = in.peek();
  if (c == '\0')
    return EOF;
  if (c == '#')
    return COMMENT;
  else if (c == '\n')
    return BLANK;
  else if (c == ' ')
    return ALIGN;
  else 
    return SEQUENCE;
}

void SRSPairParser::AssignLocations(std::string &refSeq, ssize_t refPos,
				   std::string &alignSeq, 
				   std::string &qrySeq, ssize_t qryPos,
				   ssize_t *locations, 
				   ssize_t &refGaps, ssize_t &qryGaps,
				   ssize_t refStart, ssize_t qryStart) {
  //UNUSED// char a;
  ssize_t i;
  i = 0;
  refGaps = 0;
  qryGaps = 0;
  for (; i < refSeq.size(); i++) {
    if (refSeq[i] == '-') {
      qryPos++;
      refGaps++;
    }
    else if (qrySeq[i] == '-') {
      refPos++;
      qryGaps++;
    }
    else {
      locations[refPos - refStart] = qryPos - qryStart;
      refPos++;
      qryPos++;
    }
  }
}
				   

ssize_t SRSPairParser::ParseAlignmentTriplet(std::ifstream &in, ssize_t *locations, ssize_t &refStart, ssize_t &qryStart, 
					 ssize_t &totRefGaps, ssize_t &totQryGaps) {

  std::string refStr, qryStr, alignStr, temp;
  std::string refLine, qryLine, alignLine;
  std::istringstream refStream, qryStream, alignStream;
  if (!GetLine(in, refLine)) return 0;
  if (!GetLine(in, alignLine)) return 0;
  if (!GetLine(in, qryLine)) return 0;

  refStream.str(refLine);
  qryStream.str(qryLine);
  alignStream.str(alignLine);
  
  ssize_t refPos, qryPos;
  ssize_t refGaps, qryGaps;

  refStream >> temp >> refPos >> refStr;
  alignStream >> alignStr;
  qryStream >> temp >> qryPos >> qryStr;

  if (refStart < 0) 
    refStart = refPos;

  if (qryStart < 0)
    qryStart = qryPos;

  AssignLocations(refStr, refPos,
		  alignStr, 
		  qryStr, qryPos,
		  locations, refGaps, qryGaps, refStart, qryStart);

  totRefGaps += refGaps;
  totQryGaps += qryGaps;
	return 1;
}

void SRSPairParser::FixLocations(ssize_t *&locations, ssize_t &length, ssize_t &refGaps) {
  ssize_t *newLocations = new ssize_t[length - refGaps];
  ssize_t i;
  length = length - refGaps;
  for (i = 0; i < length; i++) 
    newLocations[i] = locations[i];

  delete []locations;
  locations = newLocations;
}

ssize_t SRSPairParser::ParseSRSFile(std::string fileName, 
				 EmbossAlignment &alignment) {

  std::ifstream srsFile;

  openck(fileName, srsFile);
  ssize_t lineType;
  //UNUSED// ssize_t commentType;
  std::string line;
  //UNUSED// ssize_t num, denom;
  //UNUSED// double frac;
  // terible error handling, but at least it exists.
  if (!GetLine(srsFile, line, 16)) return 0;
  if (!GetLine(srsFile, line)) return 0;
  if (!ParseNumericLine(line, alignment.length)) return 0;
  if (!GetLine(srsFile, line)) return 0;
  if (!ParseFractionLine(line, alignment.identEq, alignment.identTot, alignment.identity)) return 0;
  if (!GetLine(srsFile, line)) return 0;
  if (!ParseFractionLine(line, alignment.simEq, alignment.simTot, alignment.similarity)) return 0;
  if (!GetLine(srsFile, line)) return 0;
  if (!ParseFractionLine(line, alignment.gapEq, alignment.gapTot, alignment.gaps)) return 0;
  if (!GetLine(srsFile, line)) return 0;
  if (!ParseNumericLine(line, alignment.alignScore)) return 0;
  if (!GetLine(srsFile, line, 4)) return 0; // Skip to the aligned portion.
  alignment.pctIdentity = double(alignment.identEq) / alignment.identTot;
  alignment.similarity  = double(alignment.simEq) / alignment.simTot;
  
  alignment.locations = new ssize_t[alignment.length];
  ssize_t i;
  for (i = 0; i < alignment.length; i++) {
    alignment.locations[i] = -1;
  }
  ssize_t refGaps = 0;
  ssize_t qryGaps = 0;
  while (srsFile) {
    if (! ParseAlignmentTriplet(srsFile, alignment.locations, alignment.refStart, alignment.qryStart, refGaps, qryGaps)) return 0;
    srsFile.get(); // get the next newline
    lineType = GetLineType(srsFile);
    switch(lineType) {
    case EOF:
      break;
    case COMMENT:
      break;
    case BLANK:
      break;
    }
  }
  FixLocations(alignment.locations, alignment.length, refGaps);
  alignment.refStart -= 1;
  alignment.qryStart -= 1;
  return 1;
}

  
