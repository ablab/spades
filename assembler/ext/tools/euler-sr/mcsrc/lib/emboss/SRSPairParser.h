/***************************************************************************
 * Title:          SRSPairParser.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef SRSPAIR_PARSER_H_
#define SRSPAIR_PARSER_H_

#include <string>
#include "EmbossAlignment.h"

class SRSPairParser {
public:
  static const ssize_t LENGTH;
  static const ssize_t IDENTITY;
  static const ssize_t SIMILARITY;
  static const ssize_t GAPS;
  static const ssize_t SCORE;

  static const ssize_t COMMENT, BLANK, SEQUENCE, ALIGN;
  static ssize_t ParseSRSFile(std::string fileName, EmbossAlignment &alignment);
  static ssize_t GetLineType(std::ifstream &in);
  static ssize_t GetCommentLineType(std::string &line);
  static ssize_t GetLine(std::ifstream &in, std::string &line, ssize_t repeat=1);
  static ssize_t ParseFractionLine(std::string &line, ssize_t &num, ssize_t &denom, double &frac);

  template <typename t>
  static ssize_t ParseNumericLine(std::string &line, t & num) {
    std::istringstream  lineStream(line);
    std::string temp;
    //UNUSED// char c;
    
    if ( ! (lineStream >> temp >> temp)) return 0;
    if (! (lineStream >> num ) ) return 0;
    
    return 1;
  }
  static ssize_t ParseAlignmentTriplet(std::ifstream &in, 
				   ssize_t *locations, 
				   ssize_t &refStart, ssize_t &qryStart,
				   ssize_t &totRefGaps, ssize_t &totQryGaps);

  static void AssignLocations(std::string &refSeq, ssize_t refPos,
			     std::string &alignSeq, 
			     std::string &qrySeq, ssize_t qryPos,
			     ssize_t *locations, 
			     ssize_t &refGaps, ssize_t &qryGaps,
			     ssize_t refStart, ssize_t qryStart);

  static void FixLocations(ssize_t *&locations, ssize_t &length, ssize_t &refGaps);
};
#endif
