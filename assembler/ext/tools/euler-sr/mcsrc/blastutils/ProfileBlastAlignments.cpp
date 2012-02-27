/***************************************************************************
 * Title:          ProfileBlastAlignments.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "parseblast/BlastParser.h"
#include "parseblast/BlastResult.h"
#include "SeqReader.h"
#include "ParseTitle.h"
#include "DNASequence.h"
#include "SeqUtils.h"
#include "utils.h"
#include <map>


typedef std::string PrefStr;

typedef std::map<std::string, ssize_t> CountMap;

ssize_t CopyPref(unsigned char* string, ssize_t *alignment, ssize_t pos, std::string &pref, ssize_t prefLen);
void IncrementPref(CountMap &map, PrefStr pref);

int main(int argc, char* argv[]) {
  std::string refFileName, qryFileName, blastFileName, outFileName;

  if (argc < 5) {
    std::cout << "usage: profileBlastAlignments refFile qryFile blastFile outFile [options]" << std::endl
							<< "      [-printErrorPrefix] Print the 4 nucleotdes before any error " << std::endl
							<< "      [-printCoords]      Print the coordinates of the top blast hit " << std::endl
							<< "      [-printStat]        Print the alignment statistics for each alignment" << std::endl
							<< "      [-printMismatch]    Print titles of queries that are not a perfect match."
							<< std::endl
							<< std::endl;
    exit(1);
  }

  refFileName = argv[1];
  qryFileName = argv[2];
  blastFileName = argv[3];
  outFileName = argv[4];
  int argi = 5;
  ssize_t printErrorPrefix = 0;
  ssize_t printCoords = 0;
  ssize_t printSummary = 1;
	ssize_t printStat = 0;
	ssize_t printMismatch = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-printErrorPrefix") == 0) {
      printErrorPrefix = 1;
    }
    if (strcmp(argv[argi], "-printCoords" ) == 0) {
      printCoords = 1;
      printSummary = 0;
    }
		if (strcmp(argv[argi], "-printStat" ) == 0) {
			printStat = 1;
		}
		if (strcmp(argv[argi], "-printMismatch") == 0 ) {
			printMismatch = 1;
		}
    ++argi;
  }
  std::ofstream summaryOut;
  openck(outFileName, summaryOut, std::ios::out);
  DNASequence refSeq, refRC;
  SeqReader::GetSeq(refFileName, refSeq, SeqReader::noConvert);
  MakeRC(refSeq, refRC);

  DNASequenceList queries;
  ReadDNASequences(qryFileName, queries);

  BlastResult blastResult;
  
  if (!(ReadBlastFile(blastFileName, blastResult))) {
    std::cout << "Could not parse the blast file " << std::endl;
  }

  ssize_t queryIndex, blastIndex;
  queryIndex = 0;
  blastIndex = 0;
  ssize_t numBadSeq = 0;
  ssize_t totMis = 0;
  ssize_t totDel = 0;
  ssize_t totIns = 0;
  ssize_t totSeq = 0;
  ssize_t totNotAligned = 0;
  std::string pref;
  pref[4] = '\0';
  CountMap misMap, insMap, delMap, indelMap;
  
  while (queryIndex < queries.size() and
	 blastIndex < blastResult.size()) {
  // The blast results must be a subset of the queries in the same order.
    std::string queryName = queries[blastIndex].namestr;
		std::string queryTitle;
		ParseTitle(queryName, queryTitle);

    if ( queryTitle != blastResult[blastIndex]->queryName) {
      queryIndex++;
      continue;
    }
    if (blastResult[blastIndex]->hsps.size() > 0) {
			//			std::cout << "q:" << queries[blastIndex].namestr << ", " << queryTitle << ", " << blastResult[blastIndex]->queryName << std::endl; 			
      ssize_t nIns, nDel, nMis;
      nIns = 0; nDel = 0; nMis = 0;
      BlastQueryMatch *query = blastResult[blastIndex];
      totSeq++; 
			ssize_t refEnd, refPos;
			refEnd = query->hsps[0].refEnd;
			refPos = query->hsps[0].refPos;

				//      if ( (refEnd > refPos  and  
				//						(double(refEnd - refPos ) / 
				//						 double(queries[queryIndex].length)  < 0.7)) or
				//					 (refPos > refEnd  and  
				//						(double(refPos - refEnd ) / 
				//						 double(queries[queryIndex].length)  < 0.7)))	{ 

			if ( abs(refEnd - refPos)
					 < 0.7 * queries[queryIndex].length ) {
				std::cout << queryTitle << " " << queries[queryIndex].length << " " << query->hsps[0].refEnd << " " << query->hsps[0].refPos << " " << query->hsps[0].refEnd - query->hsps[0].refPos << std::endl;
        totNotAligned++;
      }
      else {
      if (query->hsps[0].identity != 100) {
        // found a gapped match.  Profile the match

				ssize_t alignPos = 0;
				ssize_t qryPos = query->hsps[0].qryPos;
				unsigned char* seq, *qrySeq;
				if (query->hsps[0].strand == 0) 
					seq = refSeq.seq;
				else {
					FixMinusCoordinates (query->hsps[0], refSeq.length);
					seq = refRC.seq;
				}
		       
				qrySeq = queries[queryIndex].seq;
				nIns = nDel = nMis = 0;
				ssize_t *alignment = query->hsps[0].locations;

				/*
					Badly coded print to check out the aligned strings
					char tmpQrySeq[1000];
					char tmpRefSeq[1000];
					ssize_t ql, rl, p;
					ql = 0;
					for (p = 0; p < query->hsps[0].alignmentLength; p++ ) 
					if (alignment[p] != -1) ql++;
	
      
					strncpy(tmpQrySeq, (const char*)(qrySeq + qryPos), ql);
					rl = alignment[query->hsps[0].alignmentLength-1] - alignment[0] + 1;
					strncpy(tmpRefSeq, (const char*) &(seq[alignment[0]]), rl);
	
					tmpQrySeq[ql] ='\0';
					tmpRefSeq[rl] ='\0';
					std::cout << "aligned: " << std::endl;
					std::cout << tmpQrySeq << std::endl;
					std::cout << tmpRefSeq << std::endl;
				*/
				for (alignPos = 1; alignPos < query->hsps[0].alignmentLength; alignPos++) {
					if (alignment[alignPos] == -1) {
						++nIns;
						if (CopyPref(&qrySeq[qryPos], alignment, alignPos-1, pref, 4)) {
							IncrementPref(insMap, pref);
							IncrementPref(indelMap, pref);
						}
					}
					else {
						// If we have advanced to a deletion, first count the deletion
						// then check for mismatch
						if (alignment[alignPos] != -1 and
								alignment[alignPos-1] != -1 and
								alignment[alignPos] != alignment[alignPos-1] + 1) {
							nDel += alignment[alignPos] - alignment[alignPos-1] - 1;
							if (CopyPref(&qrySeq[qryPos],alignment, alignPos-1, pref, 4)) {
								IncrementPref(delMap, pref);
								IncrementPref(indelMap, pref);
							}
						}
						// Check for mismatch
						if (unmasked_nuc_index[qrySeq[qryPos + alignPos]] !=
								unmasked_nuc_index[seq[alignment[alignPos]]]) {
							++nMis;
							if (CopyPref(&qrySeq[qryPos], alignment, alignPos-1, pref, 4)) {
								IncrementPref(misMap, pref);
							}
						}
					}
				}
				numBadSeq++;
				if (printStat) {
					std::cout << queryTitle << " " << nIns << " " << nDel << " " << nMis << std::endl;
				}
      }
      if (printSummary) {
				summaryOut << queries[queryIndex].namestr 
									 << " "<<   query->hsps[0].identity << " " << query->hsps[0].alignmentLength 
									 << " " << nMis << " " << nIns << " " << nDel << std::endl;
      }
      else if (printCoords) {
				summaryOut << queries[queryIndex].namestr  << " "
									 << query->hsps[0].refPos << " " << query->hsps[0].refEnd << std::endl;
      }



      totMis += nMis;
      totIns += nIns;
      totDel += nDel;
    }
    }
    queryIndex++;
    blastIndex++;
  }
  std::cout << "summary: out of " << totSeq << ", " << numBadSeq 
	    << " sequences had errors and " << totNotAligned 
            << " were not aligned " << std::endl;
  std::cout << "a total of " << totMis << " mismatches " << totIns 
						<< " insertions and " << totDel << " deletions" << std::endl;
  if (printErrorPrefix) {
    std::cout << "Mutation prefix: " << misMap.size() << std::endl;
    CountMap::iterator it;
    for (it = misMap.begin(); it != misMap.end(); ++it) {
      std::cout << it->first << " " << it->second << std::endl;
    }
    std::cout << "Indel prefix: " << indelMap.size() << std::endl;
    for (it = indelMap.begin(); it != indelMap.end(); ++it) {
      std::cout << it->first << " " << it->second << std::endl;
    }
  }
  return 0;
}

ssize_t CopyPref(unsigned char* seq, ssize_t*alignment, ssize_t pos, std::string &pref, ssize_t prefLen) {
  ssize_t i;
  if(pos - prefLen < 0)
    return 0;
  
  pref = "";
  for (i = prefLen-1; i > 0; i--) {
    if (alignment[pos-i] != alignment[pos-i+1]-1)
      return 0;
    pref.push_back(seq[pos-i]);
  }
  pref.push_back(seq[pos]);
  return 1;
}

void IncrementPref(CountMap &countMap, PrefStr pref) {
  CountMap::iterator it;
  if (countMap.find(pref) == countMap.end()) {
    countMap[pref] = 0;
  }
  countMap[pref]++;
  //  std::cout << "incrementing " << pref << " " << countMap.size() << std::endl;
}
