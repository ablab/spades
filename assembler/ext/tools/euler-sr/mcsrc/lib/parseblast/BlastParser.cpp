/***************************************************************************
 * Title:          BlastParser.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/23/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BlastParser.h"
#include "BlastResult.h"
#include <fstream>
#include <assert.h>
#include "utils.h"
#include <errno.h>
#include <stdlib.h>
/*
	BlastParser::blastMajor = 0;
	BlastParser::blastMinor = 0;
	BlastParser::blastRevision = 0;
*/
void FixMinusCoordinates(BlastHSP &hsp, ssize_t seqLength) {
  ssize_t refBegin, refEnd;
  refBegin = seqLength - hsp.refPos - 1;
  refEnd   = seqLength - hsp.refEnd - 1;
  hsp.refPos = refBegin;
  hsp.refEnd = refEnd;

  ssize_t i;
  for (i = 0; i < hsp.alignmentLength; i++ ) {
    if (hsp.locations[i] != -1)
      hsp.locations[i] = seqLength - hsp.locations[i] - 1;
  }
}

ssize_t BlastParser::ParseBlastFile(std::istream &stream, BlastResult &result) {

  LineCountedIStream in(&stream);
	blastMajor = 0;
	blastMinor = 0;
	blastRev   = 0;
  std::string word, keyword;
  // start off
  //UNUSED// BlastQueryMatch *queryMatch;

  std::string kwBlastn = "BLASTN";
  std::string kwQuery  = "Query=";
  std::string kwDatabase= "Database:";
  std::string kwNoHit = "*****";
  std::string kwHits = "Sequences";
  std::string kwHSP  = "Score";
  ssize_t status = 1;
  std::string lastLine;
	keyword = "";
	ssize_t iter = 0;
  while(in.good()) {
		++iter;
		/*
		if (iter % 1000 == 0)
			std::cout << "iter; " << iter << std::endl;
		*/
		if (keyword == kwDatabase) 
		  break;
		if (keyword != kwBlastn) 
			// only read the next keyword if the previous run didn't end on the 
			// blastn keyword
			if (!(in.Stream() >> keyword)) 
				break;
    if (keyword == kwBlastn) {
      if (!ReadHeader(in)) break;
    }
    else {
      std::cout << "LINE: " << in.lineNumber << " bad format, expected " << kwBlastn << std::endl;
      in.GetLine(lastLine);
      status = 0;
      break;
    }
    BlastQueryMatch *queryMatch = new BlastQueryMatch;
    result.push_back(queryMatch);
    if (!(in.Stream() >> keyword)) 
      break;
    if (keyword == kwQuery) {
      if (!ParseQueryName(in, queryMatch->queryName)) break;

			if (queryMatch->queryName == "49") {
				Pause();
			}
    }
    else {
      std::cout << "LINE: " << in.lineNumber << " bad format, expected " << kwQuery << std::endl;
      in.GetLine(lastLine);
      status = 0;
      break;
    }

    if (!(in.Stream() >> keyword)) 
      break;
    if (keyword == kwDatabase) {
      if (!SkipDatabase(in)) break;
      
    }
    else {
      std::cout << "bad format, expected " << kwDatabase << std::endl;
      in.GetLine(lastLine);
      status = 0;
      break;
    }

    if (!(in.Stream() >>keyword)) break;
    if (keyword == kwNoHit) {
      if (!SkipEmptyHit(in)) break;
      in.Stream() >> keyword;
    }
    else {
      if (keyword != kwHSP) {
				std::cout << "LINE: " << in.lineNumber<< " bad format, expected: " << kwHSP;
				in.GetLine(lastLine);
				status = 0;
				break;
      }
			if (!SkipDatabaseHits(in)) break;
      // there were hits, parse them
			std::string sbjct;
      do {
				if (!ParseSbjct(in, sbjct)) break;
				do {
					if (in.Stream().peek() == '>' or 
							!(in.Stream() >> keyword)) {
						break;
					}
					if (keyword == kwHSP) {
						//	  std::cout << "attempting to parse a hsp!!!" << std::endl;
						BlastHSP hsp;
						if (ParseHSP(in, hsp)) {
							hsp.sbjct = sbjct;
							queryMatch->hsps.push_back(hsp);
						}
						else {
							std::cout << "could not parse hsp at " << in.lineNumber << std::endl;
							in.get();
							in.GetLine(lastLine);
							status = 0;
							break;
						}
					}
				} while (keyword == kwHSP);
			} while (in.Stream().peek() == '>');
    }
		if (blastRev <= 13) {
			if (keyword != kwDatabase) {
				std::cout << "expected postamble beginning with " << kwDatabase
									<< std::endl;
			}
			else {
				SkipPostamble(in);
			}
		}
  }
  if (status == 0) {
    std::cout << "ended parsing on line: " << lastLine << std::endl;
  }

	return 1; // success
	// TODO: return 0 on fail
}

ssize_t BlastParser::ParseAlignmentTriplet(LineCountedIStream &in,
																			 std::string &qryLine, 
																			 std::string &sbjctLine, 
																			 ssize_t &qryStart, ssize_t &qryEnd, 
																			 ssize_t &sbjctStart, ssize_t &sbjctEnd) {

  // Predcondition: the next word on input is "Query:"
  std::string keyword, word, line;
  
  if (!(in.Stream() >> keyword >> qryStart >> qryLine >> qryEnd))
    return 0;

  in.GetLine(line);
  in.GetLine(line);
  
  in.Stream() >> word >> sbjctStart >> sbjctLine >> sbjctEnd;
  if (!(in.GetLine(line))) return 0;
  if (!(in.GetLine(line))) return 0;
  if (!(in.GetLine(line))) return 0;
  return 1;
}

ssize_t BlastParser::CountQryLength(std::string &qryLine) {
  ssize_t i, length;
  char *qry;
  qry = (char*) qryLine.c_str();
  length = qryLine.length();
  ssize_t qryLength = 0;
  for (i = 0; i < length; i++ ){ 
    if (qry[i] != '-') qryLength++;
  }
  return qryLength;
}

ssize_t BlastParser::TripletToAlignment(std::string &fullQryLine, 
																		std::string &fullSbjctLine,
																		ssize_t sbjctBegin,
																		ssize_t strand,
																		ssize_t *&alignment,
																		ssize_t &qryLength) {
  assert(fullQryLine.size() == fullSbjctLine.size());

  qryLength = CountQryLength(fullQryLine);
  alignment = (ssize_t*) new ssize_t[qryLength];
  ssize_t i;
  //UNUSED// ssize_t qryStrPos = 0;
  //UNUSED// ssize_t sbjctStrPos = 0;
  ssize_t qryAlignPos = 0;
  ssize_t sbjctAlignPos = 0;
  for (i = 0; i < qryLength; i++) alignment[i] = -1;
  ssize_t alignPos = 0;
  char *qryStr, *sbjctStr;
  ssize_t alignLength;
  alignLength = fullQryLine.length();
  qryStr   = (char*) fullQryLine.c_str();
  sbjctStr = (char*) fullSbjctLine.c_str();

  while (alignPos < alignLength) {
    assert(qryAlignPos < qryLength);
    if (qryStr[alignPos] != '-' and sbjctStr[alignPos] != '-') {
      alignment[qryAlignPos] = sbjctAlignPos + sbjctBegin;
      qryAlignPos++;
      if (strand == 0)
				sbjctAlignPos++;
      else
				sbjctAlignPos--;
    }
    else if (qryStr[alignPos] == '-') {
      if (strand == 0)
				sbjctAlignPos++;
      else
				sbjctAlignPos--;
    }
    else if (sbjctStr[alignPos] == '-') {
      qryAlignPos++;
    }
    else {
      std::cout << "BLAST Parsing error: two gaps are aligned " << std::endl;
    }
    alignPos++;
  }
  assert(qryAlignPos == qryLength);
  return 1;
}
  
ssize_t BlastParser::ParseBlastVersion(std::string &line) {
	std::stringstream vstrm(line);
	std::string temp;
	// parse a version string of the format:
	// BLASTN 2.2.16 [Mar-25-2007]
	vstrm >> blastMajor;
	vstrm.get(); // get '.'
	vstrm >> blastMinor;
	vstrm.get(); // get '.'
	vstrm >> blastRev;
	return blastMajor;
}

ssize_t BlastParser::ReadHeader(LineCountedIStream &in) {
  std::string line;
  in.GetLine( line);; // BlastN 
	ParseBlastVersion(line);
  in.GetLine( line);; // blank
  in.GetLine( line);; // blank
  in.GetLine( line);; // Reference: Altschul, Stephen F., Thomas L. Madden, Alejandro A. Schaffer, 
  in.GetLine( line);; // Jinghui Zhang, Zheng Zhang, Webb Miller, and David J. Lipman (1997), 
  in.GetLine( line);; // "Gapped BLAST and PSI-BLAST: a new generation of protein database search
  in.GetLine( line);; // programs",  Nucleic Acids Res. 25:3389-3402.
  in.GetLine( line);
  //  std::cout << "read header ended on: " << line << std::endl;
	return in.good();
}

ssize_t BlastParser::ParseQueryName(LineCountedIStream &in, std::string &queryName) {
  std::string line;
  if (!(in.Stream() >> queryName)) return 0;
  in.GetLine( line);; // polish off the line;
  in.GetLine( line);; // (xx letters)
  in.GetLine( line);; // blank
  /*
    std::cout << "pqn ended on " << line << std::endl;
    std::cout << "got name: " << queryName << std::endl;
  */
  if (in.good()) {
    return 1;
  }
  else 
    return 0;
}

ssize_t BlastParser::SkipDatabase(LineCountedIStream &in) {
  std::string line;
  in.GetLine( line);; // Database: blah
  in.GetLine( line);; // 1 sequences ... blah
  in.GetLine( line);; // blank
  in.GetLine( line);; // searching done
  in.GetLine( line);
  //  std::cout << "skip database ended on " << line << std::endl;
  return 1;
}

ssize_t BlastParser::SkipEmptyHit(LineCountedIStream &in) {
  std::string line;
  in.GetLine( line);
  in.GetLine( line);
	return 1;
}

ssize_t BlastParser::SkipDatabaseHits(LineCountedIStream &in) {
  std::string line;
  in.GetLine( line);; // Score E
  in.GetLine( line);; // Sequences producing significant...
  in.GetLine( line);; // blank
  std::string match;
  do {
    in.GetLine( line);;
    //    std::cout << "got significant sequence: " << line << std::endl;
  } while (in.Stream() and in.Stream().peek() != '>');
  if (!in.good()) return 0;
  return 1;
}

ssize_t BlastParser::ParseSbjct(LineCountedIStream &in, std::string &sbjct) {
  std::string word, line;
  in.GetLine(line);
	std::stringstream sbjctLine(line);
	
  if (line[0] != '>')
    return 0;
	sbjctLine.get();
	sbjctLine >> sbjct;
  do {
    if (!(in.Stream() >> word)) {
      std::cout << "BLAST Parsing error, looking for \"Length =\"" << std::endl;
      exit(1);
    }
    in.GetLine(line);
    //    std::cout << "skipping match: " << line << std::endl;
  } while (in.good() and word != "Length");
  in.GetLine( line); // discard the next line
  //  std::cout << "ended skipmatch on: " << line << std::endl;
  // should be on a hsp
  if (!in.good()) return 0;
  return 1;
}

ssize_t BlastParser::SkipMatch(LineCountedIStream &in) {
  std::string word, line;
  in.GetLine(line);
  if (line[0] != '>')
    return 0;
  do {
    if (!(in.Stream() >> word)) {
      std::cout << "BLAST Parsing error, looking for \"Length =\"" << std::endl;
      exit(1);
    }
    in.GetLine( line);
    //    std::cout << "skipping match: " << line << std::endl;
  } while (in.good() and word != "Length");
  in.GetLine( line); // discard the next line
  //  std::cout << "ended skipmatch on: " << line << std::endl;
  // should be on a hsp
  if (!in.good()) return 0;
  return 1;
}

ssize_t BlastParser::ParseHSP(LineCountedIStream &in, BlastHSP &hsp) {
  std::string word, line;
	std::stringstream linestrm;
	//	in.GetLine(line);
	//	linestrm.str(line);
	
  if (!(in.Stream() >> word)) return 0; // =
  if (!(in.Stream() >> hsp.score)) return 0;
  if (!(in.Stream() >> word)) return 0;  // bits
	if (!(in.Stream() >> word)) return 0;  // (number)
	if (!(in.Stream() >> word)) return 0;  // Expect
	if (!(in.Stream() >> word)) return 0;  // = 
	if (!(in.Stream() >> word)) return 0;

	// char *endPtr;
	hsp.eValue = 0;
	if (word.c_str()[0] == 'e')
		word.insert(0, "1");

	if (sscanf(word.c_str(), "%Lf", &hsp.eValue) != 1) {
		std::cout << "bad conversion" << std::endl;
		return 0;
	}

	/*
	hsp.eValue = strtold(word.c_str(), &endPtr);
	if (endPtr == word.c_str()) {
		ssize_t error = errno;
		std::cout << "bad conversion with error: " << error << std::endl;
		if (error != ERANGE) 
			return 0;
	}
	*/

  in.GetLine( line); //finish off the line
  if (!in.good()) return 0;
  if (!(in.Stream() >> word // identities
				>> word // =
				>> word)) return 0;// match/leng
  in.get(); in.get(); // get 2 characters
  if (!(in.Stream() >> hsp.identity)) return 0;

  in.GetLine( line);
  /*
    std::cout << "Read hsp with score " << hsp.score << " evalue: " 
    << hsp.eValue << std::endl;
  */
  if (!(in.Stream() >> word)) return 0;
  if (!(in.Stream() >> word >> word >> word >> word)) return 0;
  if (word == "Plus") 
    hsp.strand = 0;
  else
    hsp.strand = 1;
  in.GetLine( line); // finish off the line
  in.GetLine( line); // blank
  in.GetLine( line);

  std::string fullQryLine, fullSubjectLine;
  fullQryLine = "";
  fullSubjectLine = "";
  std::string qryAlign, align, sbjctAlign;
  ssize_t tQryStart, tQryEnd, tSbjctStart, tSbjctEnd;
  ssize_t qryStart, qryEnd, sbjctStart, sbjctEnd;
  qryStart = qryEnd = sbjctStart = sbjctEnd = -1;

  while(in.Stream().peek() == 'Q') {
    ParseAlignmentTriplet(in, qryAlign, sbjctAlign, 
													tQryStart, tQryEnd, tSbjctStart, tSbjctEnd);
    /*
      std::cout << "qryAlign: '"<<qryAlign<<"'"<<std::endl;
      std::cout << "sbjAlign: '"<<sbjctAlign<<"'"<<std::endl;
    */
    if (qryStart == -1) {
      qryStart = tQryStart;
      sbjctStart = tSbjctStart;
    }
    qryEnd   = tQryEnd;
    sbjctEnd = tSbjctEnd;

    fullQryLine += qryAlign;
    fullSubjectLine += sbjctAlign;
  }
  /*
		std::cout << "got full query: " << fullQryLine << std::endl;
  */
  hsp.qryPos   = qryStart - 1;
  hsp.qryEnd   = qryEnd - 1;
  hsp.refPos   = sbjctStart - 1;
  hsp.refEnd   = sbjctEnd - 1;
  if (hsp.refEnd == -1) {
    std::cout << "error, refEnd negative " << std::endl;
  }
  /*
    std::cout << "got coordinates: " << hsp.qryPos << " " << hsp.qryEnd
    << " " << hsp.refPos << " " << hsp.refEnd << std::endl;
  */
  TripletToAlignment(fullQryLine, fullSubjectLine, 
										 hsp.refPos, hsp.strand, hsp.locations, hsp.alignmentLength);

  if (!in.Stream()) return 0;

  return 1;
}



ssize_t BlastParser::SkipPostamble(LineCountedIStream &in) {
  /* Skip:
		 Database: spneumoniae.fasta
		 Posted date:  Mar 15, 2007 10:00 AM
		 Number of letters in database: 2,160,837
		 Number of sequences in database:  1
  
		 Lambda     K      H
		 1.37    0.711     1.31 

		 Gapped
		 Lambda     K      H
		 1.37    0.711     1.31 


		 Matrix: blastn matrix:1 -3
		 Gap Penalties: Existence: 5, Extension: 2
		 Number of Hits to DB: 48
		 Number of Sequences: 1
		 Number of extensions: 48
		 Number of successful extensions: 1
		 Number of sequences better than 1.0e-06: 1
		 Number of HSP's better than  0.0 without gapping: 1
		 Number of HSP's successfully gapped in prelim test: 0
		 Number of HSP's that attempted gapping in prelim test: 0
		 Number of HSP's gapped (non-prelim): 1
		 length of query: 67
		 length of database: 2,160,837
		 effective HSP length: 13
		 effective length of query: 54
		 effective length of database: 2,160,824
		 effective search space: 116684496
		 effective search space used: 116684496
		 T: 0
		 A: 0
		 X1: 11 (21.8 bits)
		 X2: 15 (29.7 bits)
		 S1: 12 (24.3 bits)
		 S2: 24 (48.1 bits)
  */

  ssize_t i;
  std::string line;
  for (i = 0; i < 37; i++ ) {
    in.GetLine( line);
    if (!in.Stream()) return 0;
  }
  return 1;
}

ssize_t ReadBlastFile(std::string &blastFileName, BlastResult &result) {
  std::ifstream stream;
  openck(blastFileName, stream, std::ios::in);
  return BlastParser::ParseBlastFile(stream, result);
}


ssize_t ReadBlastTable(std::string &blastTableName, BlastResult &result) {
	std::ifstream stream;
  openck(blastTableName, stream, std::ios::in);
  return BlastParser::ParseBlastTable(stream, result);
}
	

ssize_t BlastParser::ParseBlastTable(std::istream &stream, BlastResult &result) {
	std::string query = "";
	std::string curQuery = "";
  LineCountedIStream in(&stream);
	ssize_t done = 0;
	while (in.Stream().good() && !done) {
		if (query == "")
			in.Stream() >> query;
		BlastQueryMatch *queryMatch;
		
		if (query != curQuery) {
			queryMatch = new BlastQueryMatch;;
			queryMatch->queryName = query;
			result.push_back(queryMatch);
			curQuery = query;
		}
		while (query == curQuery and !done) {
			BlastHSP hsp;
			//UNUSED+// ssize_t  nmatch;
			ssize_t nmismatch, ngap;
			ssize_t failed = 0;
			if (!(in.Stream() >> hsp.sbjct)) failed = 1;
			else if (!(in.Stream() >> hsp.identity)) failed = 1;
			else if (!(in.Stream() >> hsp.length)) failed = 1;
			else if (!(in.Stream() >> nmismatch)) failed = 1;
			else if (!(in.Stream() >> ngap)) failed = 1;
			else if (!(in.Stream() >> hsp.qryPos)) failed = 1;
			else if (!(in.Stream() >> hsp.qryEnd)) failed = 1;
			else if (!(in.Stream() >> hsp.refPos)) failed = 1;
			else if (!(in.Stream() >> hsp.refEnd)) failed = 1;
			else if (!(in.Stream() >> hsp.eValue)) failed = 1;
			else if (!(in.Stream() >> hsp.score)) failed = 1;

			if (failed) {
				std::cout << "Error parsing table at line: " << in.lineNumber;
				assert(0);
			}
			in.GetLine();
			queryMatch->hsps.push_back(hsp);
			if (!(in.Stream() >> query))
				done = 1;
		}
	}
	return 1; // read successfully
}
