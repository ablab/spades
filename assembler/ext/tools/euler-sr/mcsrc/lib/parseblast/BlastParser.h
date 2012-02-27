/***************************************************************************
 * Title:          BlastParser.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BLAST_PARSER_H_
#define BLAST_PARSER_H_

#include <string>
#include <iostream>

#include "BlastResult.h"
#include "BlastHSP.h"

class LineCountedIStream {
public:
  std::istream *stream;
  ssize_t lineNumber;
  LineCountedIStream(std::istream *in) {
    stream = in;
		lineNumber = 0;
  }
  std::istream &Stream() {
    return *stream;
  }
  ssize_t GetLine(std::string &line) {
    std::getline(*stream, line);
    lineNumber++;
    return (good());
  }
	char get() {
		char result;
		if ((result = stream->get()) == '\n')
			lineNumber++;
		return result;
	}
	ssize_t GetLine() {
		std::string line;
		std::getline(*stream, line);
		lineNumber++;
		return good();
	}

  ssize_t good() {
    return (stream->good());
  }
};

void FixMinusCoordinates(BlastHSP &hsp, ssize_t seqLength);

ssize_t ReadBlastFile(std::string &blastFileName, BlastResult &restlt);
ssize_t ReadBlastTable(std::string &blastTableName, BlastResult &result);

namespace BlastParser {
	static ssize_t blastMajor, blastMinor, blastRev;
	ssize_t ParseSbjct(LineCountedIStream &in, std::string &sbjct);
	ssize_t ParseBlastVersion(std::string &line);
  ssize_t ParseBlastFile(std::istream &in, BlastResult &result);
	ssize_t ParseBlastTable(std::istream &in, BlastResult &result);
  ssize_t ParseAlignmentTriplet(LineCountedIStream &in, 
			    std::string &qryLine, 
			    std::string &sbjctLine, 
			    ssize_t &qryStart, ssize_t &qryEnd, 
			    ssize_t &sbjctStart, ssize_t &sbjctEnd);

  ssize_t CountQryLength(std::string &qryLine);

  ssize_t TripletToAlignment(std::string &fullQryLine, 
			 std::string &fullSbjctLine,
			 ssize_t sbjctBegin, ssize_t strand,
			 ssize_t *&alignment, ssize_t &qryLength);

  ssize_t ReadHeader(LineCountedIStream &in);

  ssize_t ParseQueryName(LineCountedIStream &in, std::string &queryName);
  
  ssize_t SkipDatabase(LineCountedIStream &in);
  ssize_t SkipEmptyHit(LineCountedIStream &in);
  ssize_t SkipDatabaseHits(LineCountedIStream &in);
  ssize_t SkipMatch(LineCountedIStream &in);
  ssize_t ParseHSP(LineCountedIStream &in, BlastHSP &hsp);
  ssize_t SkipPostamble(LineCountedIStream &in);
	ssize_t ParseBlastTable(std::istream &stream, BlastResult &result);
};


#endif
