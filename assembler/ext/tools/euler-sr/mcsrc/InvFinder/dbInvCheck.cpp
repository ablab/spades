/***************************************************************************
 * Title:          dbInvCheck.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std includes
#include <ios>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <unistd.h>

// 3rd party includes
#include "mysql/mysql.h"

// my stuff
#include "utils.h"
#include "lav/LAVUtils.h"
#include "lav/LAVReader.h"
#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVBlock.h"
#include "lav/LAVAlignedContig.h"
#include "net/NetDB.h"
#include "blocks/BlockDB.h"
#include "blocks/dbOrthoPos.h"
#include "blocks/BlockLib.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "alignutils.h"

class DBDescription {
public:
  std::string dbName;
  std::map<std::string, std::string> specTable;
};

typedef std::map<std::string, DBDescription*> DBDescriptionMap;

const ssize_t FORWARD = 0;
const ssize_t REVERSE = 1;

void PrintUsage() {
  std::cout << "dbinvcheck: a database based inversion checker "  << std::endl
	    << "(filter palindromes and off-diagonal alignments)." << std::endl;
  std::cout << "usage: dbinvcheck  [-R refFile -Q qryFile -f defaults-file -d database -s seqNameTable " << std::endl
	    << "                    -w window]  lavFile refFile qryFile chrom outfile." << std::endl;
}

void GetOverlappingCoordinates(DBQuery &query, std::string &netTableName,
			       ssize_t refLength, ssize_t qryLength,
			       ssize_t &refStart, ssize_t &refEnd, 
			       ssize_t &qryStart, ssize_t &qryEnd);


ssize_t InitEnv(int argc, char* argv[], 
	    std::string &defaultsFile,
	    std::string &dbName,
	    std::string &lavFile,
            std::string &refSpec, std::string &qrySpec, std::string &chrom,
            std::string &refFileName, std::string &qryFileName,
	    std::string &outFileName,
 	    std::string &seqTableName,
	    std::string &scoreMatFile,
	    ssize_t &explain,
	    ssize_t &window,
	    ssize_t &verbose) {
  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "d:s:w:m:vf:R:Q:x")) != EOF){
    switch(copt) {
    case 'x':
      explain = 1;
      continue;
    case 'R':
      refFileName = optarg;
      continue;
    case 'Q':
      qryFileName = optarg;
      continue;
    case 'f':
      defaultsFile = optarg;
      continue;
    case 'v':
      verbose++;
      continue;
    case 'm':
      scoreMatFile = optarg;
      continue;
    case 'd':
      dbName = optarg;
      continue;
    case 's':
      seqTableName = optarg;
      continue;
    case 'w':
      window = atoi(optarg);
      continue;
    default:
      PrintUsage();
      exit(0);
      continue;
    }
  }
  i = optind;
  if (i < argc) { 
    lavFile = argv[i++];
  }
  if (i < argc) {
    refSpec = argv[i++];
  } 
  if (i < argc) {
     qrySpec = argv[i++];
  }
  if (i < argc) {
     chrom = argv[i++];
  }
  if (i < argc) {
    outFileName = argv[i++];
  }
  else {
    PrintUsage();
    exit(1);
  }
}

void CreateTableName(std::string &refSpecies,
		     std::string &file,
		     std::string &qrySpecies,
		     std::string &tableName) {
  tableName = refSpecies + file + "_chain_" + qrySpecies;
}

void JoinQueryResult(MYSQL_RES *res,
		     ssize_t &qStart,
		     ssize_t &qEnd) {

  // Given a query, find the beginning and the ending of the query. 
  // 
  qStart = -1;
  qEnd   = -1;

  // For now assume that the result is sorted in order of query, and
  // simply use the result of the first and last row.
  ssize_t numRows = mysql_num_rows(res);
  MYSQL_ROW row;
  if (numRows > 0) {
    LAVRow lavRow;
    row = mysql_fetch_row(res);
    GetBlock(row, lavRow);
    qStart = lavRow.qStart;
    mysql_data_seek(res, numRows-1);
    row = mysql_fetch_row(res);
    GetBlock(row, lavRow);
    qEnd = lavRow.qEnd;
  }
}

ssize_t GetInversionRows(DBQuery &query, std::string tableName, 
		     ssize_t refStart, ssize_t refEnd, 
		     ssize_t qryStart, ssize_t qryEnd,
		     std::vector<LAVRow> & inversion);


int main(int argc, char* argv[]) {
  std::string posFileName;
  std::string outputFileName;
  std::string seqTableName;
  std::string dbName;
  std::string lavFileName;
  std::string scoreMatName;
  std::string refSpec, qrySpec, chrom;
  std::string defaultsFile;
  std::string tFileName, qFileName;
  LAVRow lavBlock;
  dbName       = "lav_alignments";
  seqTableName = "sequences";
  char *home = getenv("HOME");
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";

  ssize_t window = 100;
  ssize_t verbose = 0;
  defaultsFile = "";
  ssize_t explain = 0;
  InitEnv(argc, argv, 
	  defaultsFile,
	  dbName, lavFileName, refSpec, qrySpec, chrom,
          tFileName, qFileName,
          outputFileName, seqTableName,
	  scoreMatName, 
	  explain,
	  window,
	  verbose);

  FloatMatrix scoreMat;

  std::ofstream outFile;
  openck(outputFileName, outFile, std::ios::app);
  outFile << refSpec << " " << qrySpec << " " << chrom << std::endl;

  ReadScoreMatFile(scoreMatName, scoreMat);

  
  // calcluate the average penalty:
  double avgPenalty;
  ssize_t smi, smj;
  for (smi = 0; smi < 4; smi++ ) 
    for (smj = 0; smj < 4; smj++) 
      if (smi != smj) 
	avgPenalty += scoreMat[smi][smj] - scoreMat[smi][smi];

  avgPenalty = avgPenalty / 12;
  std::cout << "avg penalty: " << avgPenalty << std::endl;
  //ParseFileName(lavFileName, refSpec, qrySpec, chrom);
  
  std::string tableName, 
    chainTableName, tmpTableName, netTableName,
    refSelfTable, qrySelfTable,
    refSelfChain, qrySelfChain;
  
  refSelfTable = refSpec + "_" + chrom + "_self";
  qrySelfTable = qrySpec + "_" + chrom + "_self";

  refSelfChain = refSelfTable + "_chain";
  qrySelfChain = qrySelfTable + "_chain";

  tableName      = refSpec + "_" + qrySpec + "_" + chrom;
  chainTableName = tableName + "_chain";
  netTableName   = tableName + "_net";
  tmpTableName   = tableName + "_temp";

  if (tFileName == "")
    tFileName = refSpec + "." + chrom + ".fa";
  if (qFileName == "")
    qFileName = qrySpec + "." + chrom + ".fa";

  LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

  DNASequence tSeq, qSeq;
  DNASequence tFrag, qFrag, qrFrag;

  SeqReader::GetSeq(tFileName, tSeq, SeqReader::noConvert);
  SeqReader::GetSeq(qFileName, qSeq, SeqReader::noConvert);

  ssize_t orientation = 0;
  ssize_t reverse = 1;

  if (DetermineOrientation(&lavFile)) {
    // orientation is -, reverse is '+'
    reverse = 0;      
    orientation = 1;
  }

  // Find blocks in the reverse that do not have corresponding
  // self alignments.
  ssize_t ac, a, b;
  LAVAlignedContig *alignedContig;
  LAVBlock *block, *block2;
  ssize_t lengthAligned, numMatches, numMismatches;
  ssize_t invLengthAligned, invNumMatches, invNumMismatches;
  ssize_t refPalLengthAligned, refPalNumMatches, refPalNumMismatches;
  ssize_t qryPalLengthAligned, qryPalNumMatches, qryPalNumMismatches;
  ssize_t foundInSelf;
  // Connect to the database;
  MYSQL *dbConn;
  
  ConnectToDB(dbConn, dbName, defaultsFile);
  
  // Create the query and associate it with the databae
  DBQuery query(dbConn);
  std::vector<LAVBlock*> trueInversions;
  ssize_t refStart, refEnd, qryStart, qryEnd;
  MYSQL_RES *res;
  MYSQL_ROW row;
  ssize_t numRows;
  ssize_t qryLength, refLength;
  
  ssize_t qryForwardStart, qryForwardEnd, refReverseStart, refReverseEnd;
  // find out what chain the top level net corresponds to.


  std::vector<ssize_t> firstLevel;
  std::string firstLevelP;
  StringSet  tempTables;
  
  if (verbose == 2) {
    std::cout << "creating " << tmpTableName << std::endl;
    std::cout << query.str() << std::endl;
  }

  //  query.execute();
  // Get the length of the query string
  query << "SELECT length from " << seqTableName << " WHERE name="
	<< dbstr(refSpec) << " AND sequence=" << dbstr(chrom) ;
  GetOneInt(query, refLength);

  query << "SELECT length from " << seqTableName << " WHERE name="
	<< dbstr(qrySpec) << " AND sequence=" << dbstr(chrom) ;
  GetOneInt(query, qryLength);
  
  ssize_t refOvpStart, refOvpEnd, qryOvpStart, qryOvpEnd;

  GetOverlappingCoordinates(query, netTableName,
			    refLength, qryLength,
			    refOvpStart,  refOvpEnd, qryOvpStart,  qryOvpEnd);

  std::cout << "got overlapping coordinates: " << refOvpStart << " " 
	    << refOvpEnd << " " << qryOvpStart << " " << qryOvpEnd << std::endl;

  double invScore;
  if (verbose > 1) {
    std::cout << "processing " << lavFile.alignments.size() << " alignments " << std::endl;
  }
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    if (alignedContig->qryContig.strand == reverse) {
      foundInSelf = 0;

      for (a = 0; a < alignedContig->alignments.size(); a++ ) {
	// Try to join gapped alignments into one.
	ssize_t a2;
	block = alignedContig->alignments[a];
	refStart = block->refBegin;
	refEnd   = block->refEnd;
	qryStart = block->qryBegin;
	qryEnd   = block->qryEnd;

	if (verbose) {
	  std::cout << "starting out on " << refStart << " " << refEnd << " " 
		    << qryStart << " " << qryEnd << std::endl;
	}

	// Join all the alignments in the same net into one.

	// Step 1.  Find the net corresponding to this alignment.
	ssize_t blockId;
	query << "SELECT id from " << tableName << " where "
	      << " tStart = " <<  block->refBegin << " AND "
	      << " qStart = " <<  block->qryBegin;
	GetOneInt(query, blockId);

	// Step 2.  Find the chain corresponding to this block.
	ssize_t chainId;
	query << "SELECT chainid from " << chainTableName 
	      << " where lavid = " << blockId;
	if (! GetUpToInt(query, chainId)) 
	  chainId = -1;

	// Step 3.  Find the net that corresponds to this chain.
	ssize_t netId;
	query << "SELECT * from " << netTableName << " WHERE chainid = " << chainId;
	query.execute();
	res = mysql_store_result(query.conn);
	if (mysql_num_rows(res) <= 0) {
	  // This inversion is not part of a net.  Try to manually glue together the 
	  // pieces of the alignmetn
	  if (verbose) {
	    std::cout << " trying to join " << block->refBegin << " " << block->refEnd
		      << ", " << block->qryBegin << " " << block->qryEnd << std::endl;
	  }
	  for (a2 = a+1; a2 < alignedContig->alignments.size(); a2++) {
	    block2 = alignedContig->alignments[a2];
	    // If the blocks are adjacent
	    if ( (block2->qryBegin >= qryEnd  or  
		  (abs(qryEnd - block2->qryBegin) < 0.5 * 
		   (block2->qryEnd - block2->qryBegin))) and
		 (block2->refBegin >= (refEnd - 20) or 
		  (abs(refEnd - block2->refBegin) < 0.5 * 
		   (block2->refEnd - block2->refBegin))) and
		 ( (block2->qryBegin - qryEnd) < 1000 ) and
		 ( (block2->refBegin - refEnd) < 1000 ) ) {
	      a = a2 ; // skip considering this inversion.
	      if ( verbose ) {
		std::cout << "joining " << refStart << " " << block2->refBegin
			  << ", " << qryStart << " " << block2->qryBegin << std::endl;
	      }
	      qryEnd = block2->qryEnd;
	      refEnd = block2->refEnd;
	    }
	  }
	}
	else {
	  // The inversion is stored as a net.  Find the end of the net
	  row = mysql_fetch_row(res);
	  ssize_t netTStart, netTSize, netTEnd, netQStart, netQEnd, netQSize;
	  netTStart = atoi(row[6]);
	  netTSize  = atoi(row[7]);
	  netTEnd   = netTStart + netTSize;
	  netQStart = atoi(row[8]);
	  netQSize  = atoi(row[9]);
	  netQEnd   = netQStart + netQSize;

	  refEnd  = std::max(netTEnd, block->refEnd);
	  qryEnd  = std::max(qryLength - netQStart, block->qryEnd);
	  if (verbose) {
	    std::cout << "found inversion on net, computed qryend from " 
		      << qryLength << " - " << netQStart << " = " << qryLength - netQStart
		      << " and " << block->qryEnd << std::endl;
	  }
	  for (a2 = a+1; a2 < alignedContig->alignments.size(); a2++) {
	    block2 = alignedContig->alignments[a2];
	    if (verbose > 1) {
	      std::cout << "checking for same inv: " << block2->refEnd << " " 
			<< netTEnd << std::endl;
	    }
	    if (block2->refEnd <= netTEnd) {
	      if (verbose) {
		std::cout << "advancing from " << a << " to " << a2 << std::endl;
	      }
	      a = a2;
	      //	      break;
	    }
	  }
	}

	// The coordinates of the ref self-alignment 
	// are always the same (regardless of the orientation of the
	// qry). 
      
	if (verbose) {
	  std::cout << std::endl << std::endl << std::endl;
	  std::cout << "considering: " << refStart << " " << refEnd << " " 
		    << qryStart << "(" << qryLength - qryStart << ") " 
		    << qryEnd << "(" << qryLength - qryEnd <<  ") ..... " << std::endl;
	}

	qryForwardStart = qryLength - qryEnd;
	qryForwardEnd = qryLength - qryStart;

	refReverseStart = refLength - refEnd;
	refReverseEnd   = refLength - refStart;

	if (refStart > refOvpEnd or
	    refEnd  < refOvpStart or
	    qryForwardStart > qryOvpEnd or
	    qryForwardEnd < qryOvpStart) {
	  if (verbose) {
	    std::cout << "skipping because coordinates are outside of orthologous region " 
		      << std::endl;
	    std::cout << refStart << " " << refEnd << " " << qryForwardStart 
		      << " " << qryForwardEnd << std::endl;
	    std::cout << refOvpStart << " " << refOvpEnd << " " << qryOvpStart 
		      << " " << qryOvpEnd << std::endl;
	  }
	  if (explain) {
	    outFile << "EX: " << lavFileName << " " << refStart << " " << refEnd 
		    << " out side orthologous region " << std::endl;
	  }
	  continue;
	}
	    
	std::vector<LAVRow> inversionRows;
	GetInversionRows(query, tableName, refStart, refEnd, qryStart, qryEnd, inversionRows);
	invScore = ScoreRows(tSeq, qSeq, inversionRows, -1, -1, 
			     invLengthAligned, 
			     invNumMatches, invNumMismatches, scoreMat);

	// Get all of the rows corresponding to the inversion.

	// Check 1.  
	// Look to see if there is an alignment in the reference
	// sequence that covers the same area as the inversion. If
	// this is the case, the inversion is just a palindrome. 

	// The tmprefselftable only has the reverse strand alignments

	// It's possible that many chained alignments are present at the 
	// position of the inversion. 

#ifdef CHECK_PALINDROME_CODE
	query << "SELECT distinct(chainid) FROM " 
	      << refSelfTable << " WHERE " 
	      << " ((tStart >= " << refStart 
	      << "  AND tStart <= " << refEnd << ") "
	      << " OR (tEnd >= " << refStart
	      << "     AND tEnd <= " << refEnd << ") "
	      << " OR (tStart < " << refStart 
	      << "     AND tEnd > " << refEnd << "))";

	query.execute();
	
	res = mysql_store_result(query.conn);
	ssize_t refPalindrome = 0;
	numRows = mysql_num_rows(res);
	if (verbose > 1) {
	  std::cout << "got " << numRows << " reference palindromic rows with " 
		    << query.prevQuery << std::endl;
	}
	if (numRows > 0) {
	  while ((row = mysql_fetch_row(res)) != NULL) {
	    chainId = atoi(row[0]);
	    
	    query << "SELECT * FROM "
		  << refSelfTable << " WHERE chainid = " << chainId;
	    if (verbose == 2) {
	      std::cout << "fetching palindrome rows with: " << query.str() << std::endl;
	    }
	    query.execute();
	    MYSQL_RES *refPalRes;
	    if ((refPalRes = mysql_store_result(query.conn)) == NULL)
	      HandleError(query.conn);
	    
	    numRows = mysql_num_rows(refPalRes);
	    if (numRows > 0) {
	      std::vector<LAVRow> palindromeRows;
	      LAVRow lavRow;
	      double palindromeScore;
	      while (row = mysql_fetch_row(refPalRes)) {
		FetchRow(row, lavRow);
		palindromeRows.push_back(lavRow);
	      }
	      //	  row = mysql_fetch_row(res);
	      palindromeScore = ScoreRows(tSeq, tSeq, palindromeRows, refStart, refEnd, 
					  refPalLengthAligned, refPalNumMatches, refPalNumMismatches,
					  scoreMat);
	      if (verbose == 2) {
		std::cout << "got ref palindrome score: " << palindromeScore << " inversion score: " 
			  << invScore<< std::endl;
		std::cout << "align stats: length aligned: " << lengthAligned << " matches: " 
			  << numMatches << " mismatches " << numMismatches 
			  << " compare to " << invNumMatches << " " << invNumMismatches << std::endl;
	      }
	      // Don't bother adding this to the queue.
	      if (palindromeScore < invScore) {
		refPalindrome = 1;
		break;
	      }
	    }
	  }
	}
	if (refPalindrome) {
	  if (verbose) {
	    std::cout << "palindromic sequence in ref " << std::endl;
	  }
	  //	  continue;
	}

	// Check 2.
	// Look to see if there is an alignment in the query sequence
	// that covers the same area as that of the inversion.  If so,
	// the inversion is just a palindrome.  Most of the time this 
	// shouldn't do anything, since a palindrome in the ref strand
	// should be preserved here.
	//

	// The coordinate of the qry Self alignment are different 
	// if the orientation of the query sequence is in reverse 
	// direction from the reference.  In this case inversions are
	// on the + alignments (0 in lav format), and need to be
	// translated to starting at the 5' end (length - pos).
	// Otherwise the position does not need to be changed.


	ssize_t qryPalindrome;
	// Find all alignments that pass through this
	if (orientation == FORWARD) {
	  query << "SELECT distinct(chainid) FROM "
		<< qrySelfTable << " WHERE qStart >= " << qryStart 
		<< " AND qStart <= " << qryEnd
		<< " AND tStart >= " << qryForwardStart 
		<< " AND tEnd <= " << qryForwardEnd;

	}

	query.execute();

	if ((res = mysql_store_result(query.conn)) == NULL)
	  HandleError(query.conn);
	
	numRows = mysql_num_rows(res);
	if (verbose > 1) {
	  std::cout << "got " << numRows << " query palindromic rows with " 
		    << query.prevQuery << std::endl;
	}

	if (numRows > 0) {

	  // At least one chained inversion is present in the self copy.  Iterate through all 
	  // of them to see if they are high scoring palindromes.
	  qryPalindrome = 0;
	  while (row = mysql_fetch_row(res)) {
	    chainId = atoi(row[0]);
	    query << "SELECT * FROM " << qrySelfTable << " WHERE chainId = " << chainId
		  << " ORDER BY tStart ";
	    query.execute();

	    MYSQL_RES *qryPalChainRes;
	    if ((qryPalChainRes = mysql_store_result(query.conn)) == NULL)
	      HandleError(query.conn);

	    if (mysql_num_rows(qryPalChainRes) > 0) {
	      std::vector<LAVRow> palindromeRows;
	      double palindromeScore;
	      LAVRow lavRow;
	      while ((row = mysql_fetch_row(qryPalChainRes)) != NULL) {
		GetBlock(row, lavRow);
		palindromeRows.push_back(lavRow);
		/*
		std::cout << "qry palindrome row: " << lavRow.tStart << " " << lavRow.tEnd << std::endl;
		*/
	      }
	      if (verbose == 2) {
		std::cout << "got " << palindromeRows.size() << " query palindrome rows with " << std::endl;
		std::cout << query.prevQuery << std::endl;
	      }
	      palindromeScore = ScoreRows(qSeq, qSeq, palindromeRows, qryForwardStart, qryForwardEnd, 
					  qryPalLengthAligned, qryPalNumMatches, qryPalNumMismatches,
					  scoreMat);
	      /*
	      std::cout << "qry palindrome score: " << palindromeScore 
	      << " inversion score: " << invScore << std::endl;
	      
	      std::cout << "qry palindrom from " << palindromeRows[0].tStart << " to " 
			<< palindromeRows[palindromeRows.size() - 1].tEnd << " ... " 
			<< palindromeRows[0].qStart << " " 
			<< palindromeRows[palindromeRows.size() - 1].qEnd << std::endl;
	      std::cout << "got qry palindrome score: " << palindromeScore 
	      << " inversion score: " << invScore<< std::endl;
	      */
	      if (verbose == 2) {
		std::cout << "palindrome: " << lavBlock.tStart << " " << lavBlock.tEnd
			  << " " << lavBlock.qStart << " " << lavBlock.qEnd << std::endl;
		std::cout << "ps: " << palindromeScore << " is: " << invScore << std::endl;
		std::cout << "align stats: length aligned: " << lengthAligned << " matches: " 
			  << numMatches << " mismatches " << numMismatches << std::endl;
	      }
	      if (palindromeScore < invScore) {
		qryPalindrome = 1;
		break;
	      }
	    }
	  } // End going through chained rows.
	  // continue;
	} // End checking to see if there are any chains through the palindrome.
      
	if (qryPalindrome == 1) {
	  if (verbose) {
	    std::cout << "palindromic sequence found in qry " 
		      <<  std::endl;
	  }
	  //continue;
	}

	mysql_free_result(res);
#endif
      
	// Now look to see if any alignments happen in the same row (appear twice 
	// in the query sequence).
	// Assume that the query alignment in questin is in the reverse strand (although
	// this isn't the case for sequenes that are already in oppposite direction).
      
	// Look to see if this inversion is on the diagonal.

	// Get the beginning of the alignment that is closest to the 
	// middle of the inversion.
	LAVRow prevBlock, nextBlock;
	ssize_t maxTStart;

	prevBlock.chainId = -1;
	nextBlock.chainId = -1;
	nextBlock.tStart = nextBlock.tEnd = nextBlock.qStart = nextBlock.qEnd = -1;
	prevBlock.tStart = prevBlock.tEnd = prevBlock.qStart = prevBlock.qEnd = -1;
        if (verbose == 2) {
	   std::cout << "fetching surrounding blocks from " << tmpTableName << " " 
                     << refStart << " " << refEnd << std::endl;
        }
	FetchSurroundingBlocks(query, tmpTableName, refStart, refEnd, -1, -1, 
			       nextBlock, prevBlock, verbose);

 
	// Look to see if this inversion is about on the diagonal or
	// not.
	// If the two points of the surrounding alingments are closer
	// to eachother than the they are to the corrected inversion,
	// the inverison is off diagonal.

	double blockDist, invDistE, invDistB;

	double invSlope, refSlope, invInt, refInt;
	double qryIsect, refIsect;

	ssize_t y1, y2;
	if (verbose) {
	  std::cout << "prevblock: " << prevBlock.tEnd << ", " << prevBlock.qEnd
		    << " nextblock: " << nextBlock.tStart << ", " << nextBlock.qStart
		    << std::endl;
	}
	if (prevBlock.tEnd == -1 or
	    prevBlock.qEnd == -1 or
	    nextBlock.tStart == -1 or
	    nextBlock.tEnd == -1 or 
	    prevBlock.tEnd < refOvpStart or
	    prevBlock.qEnd < qryOvpStart or 
	    nextBlock.tStart > refOvpEnd or 
	    nextBlock.qStart > qryOvpEnd) {
	  // There are no surrounding blocks here. Keep going.
	  if (verbose) {
	    std::cout << "skipping because a surrounding alignment wasn't found "<< std::endl;	
	    if (verbose == 2) {
	      std::cout << "used query " << query.prevQuery << std::endl;
	    }
	  }
	  if (explain) {
	    outFile << "EX: " << lavFileName << " " << refStart << " " << refEnd 
		      << " surrounding alignment not found " << std::endl;
	  }

	  continue;
	}
	    
	if (qryEnd - qryStart != 0 and 
	    nextBlock.tStart - prevBlock.tEnd != 0) {
	  // Compute the equation of the line of the alignments
	  // surrounding the inversion. 
	  if (orientation == FORWARD) {
	    y1 = prevBlock.qEnd;
	    y2 = nextBlock.qStart;
	  }
	  else {
	    y1 = qryLength - prevBlock.qEnd;
	    y2 = qryLength - nextBlock.qStart;
	  }
	  if (verbose) {
	    std::cout << y2 << " " << y1 << " " << y2 - y1 << " " << nextBlock.tStart << " " 
		      << prevBlock.tEnd << " " << nextBlock.tStart - prevBlock.tEnd 
		      << std::endl;
	  }
	  refSlope  = double(y2 - y1) / (nextBlock.tStart - prevBlock.tEnd);
	  refInt = -refSlope*prevBlock.tEnd + y1;
	  
	  // Compute the equation of the inversion.
	  if (orientation == FORWARD) {
	    y1 = qryForwardEnd;
	    y2 = qryForwardStart;
	  }
	  else {
	    y1 = qryStart;
	    y2 = qryEnd;
	  }
	  if (verbose) {
	    std::cout << y2 << " " << y1 << " " << y2 - y1 << " " << refEnd << " " 
		      << refStart << " " << refEnd - refStart  
		      << std::endl;
	  }
	  invSlope = double(y2 - y1) / (refEnd - refStart);
	  invInt = -invSlope*refStart + y1;

	  refIsect = (invInt - refInt) / (refSlope -  invSlope);
	  qryIsect = refSlope * refIsect + refInt;
	  if (verbose) {
	    std::cout << "refeqn: " << refSlope << " " << refInt << std::endl;
	    std::cout << "inveqn: " << invSlope << " " << invInt << std::endl;
	    std::cout.precision(7);
	    std::cout << " intersection: " << refIsect << " " << qryIsect << std::endl;
	  }

	  double distR, distI; 
	  // distR: distance to the reverse strand from the intersection
	  // distI: distance to the inversion from the intersection
	  if (orientation == FORWARD) {
	    distR = 0.5 * 
	      sqrt((refIsect - prevBlock.tEnd) * double(refIsect - prevBlock.tEnd) + 
		   (qryIsect - prevBlock.qEnd) * double(qryIsect - prevBlock.qEnd)) +
	      sqrt((refIsect - nextBlock.tStart) * double(refIsect - nextBlock.tStart) + 
		   (qryIsect - nextBlock.qStart) * double(qryIsect - nextBlock.qStart));
	  
	    distI = sqrt((refIsect - ((refEnd + refStart)/2.0))*(refIsect - ((refEnd + refStart)/2.0))
			 + (qryIsect - ((qryForwardEnd + qryForwardStart)/2.0)) * 
			 (qryIsect - ((qryForwardEnd + qryForwardStart)/2.0)));
	  }
	  else {
	    distR = sqrt(double(refIsect - prevBlock.tEnd)
									 * double(refIsect - prevBlock.tEnd)
									 + 
									 double(qryIsect - (qryLength - prevBlock.qEnd))
									 * double(qryIsect - (qryLength - prevBlock.qEnd)));
	  
	    distI = sqrt((refIsect - ((refEnd + refStart)/2.0))
									 * (refIsect - ((refEnd + refStart)/2.0))
									 +
									 (qryIsect - (qryLength - (qryForwardEnd + qryForwardStart)/2.0))
									 * (qryIsect - (qryLength - (qryForwardEnd + qryForwardStart)/2.0)));
	  }
	  if (verbose) {
	    std::cout << "check for along diagonal: " << std::endl;
	    std::cout << "refIsect: " << refIsect << " refstart " << refStart << " refEnd " << refEnd << std::endl;
	    std::cout << "qryIsect: " << qryIsect << " qryStart " << qryForwardStart 
		      << " qryEnd " << qryForwardEnd << std::endl;
	    std::cout << " distance from intersection to middle of inversion: " << distI << std::endl;
	    std::cout << " distance from intersection to middle of reference : " << distR << " ratio: " << distI / distR 
		      << std::endl;
	  }
	  if (((refIsect >= refStart and refIsect <= refEnd) and 
	       ( qryIsect >= qryForwardStart and qryIsect <= qryForwardEnd))
	      or ( (refIsect >= prevBlock.tEnd and refIsect <= nextBlock.tStart ) and
		   (distI / distR < 1) ) ) {
	    if (verbose) {
	      std::cout << "intersection is along diagonal " << std::endl;
	    }
	  }
	  else {
	    if (verbose) {
	      std::cout << "intersection is not on the diagonal " << std::endl;
	      std::cout << "distI: " << distI << " distR: " << distR
			<< " ratio: " << distI / distR << " " << qryForwardStart << " "
			<< qryForwardEnd << std::endl;
	    }
	    if (explain) {
	      outFile << "EX: " << lavFileName << " " << refStart << " " << refEnd 
		      << " off-diagonal inversion " << std::endl;
	    }
	    continue;
	  }
	}
	  

	// Now check to see if there is a sequence corresponding to a chain that 
	// covers the position of the inversion.

	query << "SELECT " << tmpTableName << ".id, tStart, tEnd, qStart, qEnd, strand, chainId" 
	      << " from " << tmpTableName 
	      << " WHERE (( tStart >= " << refStart 
	      << " AND tStart <= " << refEnd << ")"
	      << "        OR ( tEnd >= " << refStart 
	      << " AND tEnd <= " << refEnd << ") "
	      << "        OR ( tStart <= " << refStart 
	      << " AND tEnd >= "  << refEnd << "))"
	      << " AND strand =  " << orientation << std::endl; // want to make sure this isn't a high 

	// scoring inversion.
	query.execute();
	res = mysql_store_result(query.conn);
	std::vector<LAVRow> chainRows;
	LAVRow *lavChainRow;
	LAVRow lavRow;
	while (row = mysql_fetch_row(res)) {
	  FetchRow(row, lavRow);
	  chainRows.push_back(lavRow);
	  /*
	    std::cout << "got row to compare: " << lavRow.tStart << " " << lavRow.tEnd << " " 
	    << lavRow.qStart << " " << lavRow.qEnd << std::endl;
	  */
	}
	double  chainScore;
	if (chainRows.size() > 0) {
	  // Compare the score of the potential inversion with that of the chained seq.
	  chainScore = ScoreRows(tSeq, qSeq, chainRows, refStart, refEnd, 
				 lengthAligned, numMatches, numMismatches, scoreMat);
	  if (verbose > 1) {
	    std::cout << "ref score on the chain, pos:" << refStart << " to: " << refEnd
		      << " to query : " << chainScore << " compare to inv: " 
		      << invScore << std::endl;
	    std::cout << "align stats: length aligned: " << lengthAligned << " matches: " 
		      << numMatches << " mismatches " << numMismatches 
		      << " compare to " << invLengthAligned << " " 
		      << invNumMatches << " " << invNumMismatches << std::endl;
	    std::cout << "scores: " << invScore << " " << chainScore
		      << " allowing for " << (avgPenalty*0.5*(numMismatches)) << std::endl;
	  }
	  if (invScore >= (chainScore - (avgPenalty*0.5*(numMismatches))) and
	      lengthAligned*2 > invLengthAligned) {
	  //	  if (invScore >= chainScore) {
	    if (verbose) {
	      std::cout << "better score on the chain to ref " << invScore << " " << chainScore 
			<< " " << chainScore 
			<< " allowing for " << (avgPenalty*0.5*(numMismatches)) << std::endl;
	    }
	    if (explain) {
	      outFile << "EX: " << lavFileName << " " << refStart << " " << refEnd 
		      << " better alignment on reference " << std::endl;
	    }

	    continue;
	  }
	}
	else {
	  if (verbose) {
	    std::cout << "didin't find any foward strand alignment with: "
		      << query.prevQuery << std::endl;
	  }
	}
	query << "SELECT " << tmpTableName << ".id, tStart, tEnd, qStart, qEnd, strand, chainId" 
	      << " from " << tmpTableName 
	      << " WHERE (( qStart >= " << qryForwardStart 
	      << " AND qStart <= "  << qryForwardEnd << ")"
	      << "      OR ( qEnd >= " << qryForwardStart 
	      << " AND qEnd <= " << qryForwardEnd << ")"
	      << "      OR ( qStart <= " << qryForwardStart 
	      << " AND qEnd >= " << qryForwardEnd << "))"
	      << " AND strand =  " << orientation << std::endl; // want to make sure this isn't a high 

	// scoring inversion.
	query.execute();
	res = mysql_store_result(query.conn);
	chainRows.clear();
	ssize_t qryMatchRefStart, qryMatchRefEnd;
	qryMatchRefStart = -1;
	qryMatchRefEnd   = -1;
	while (row = mysql_fetch_row(res)) {
	  FetchRow(row, lavRow);
	  if (verbose > 1) {
	    std::cout << lavRow.qStart << "  <? " << qryForwardStart << " <? " << lavRow.qEnd << std::endl;
	    std::cout << lavRow.qStart << "  <? " << qryForwardEnd << " <? " << lavRow.qEnd << std::endl;
	  }
	  if (lavRow.qStart <= qryForwardStart and lavRow.qEnd >= qryForwardStart) {
	    qryMatchRefStart = lavRow.tStart + (qryForwardStart - lavRow.qStart);
	    if (verbose > 1) {
	      std::cout << "qmrs: " << qryMatchRefStart << std::endl;
	    }
	  }
	  if (lavRow.qStart <= qryForwardEnd and lavRow.qEnd >= qryForwardEnd) {
	    qryMatchRefEnd   = lavRow.tEnd - (lavRow.qEnd - qryForwardEnd);
	    if (verbose > 1) {
	      std::cout << "qmre: " << qryMatchRefEnd << std::endl;
	    }
	  }
	  chainRows.push_back(lavRow);
	}
	if (chainRows.size() > 0) {
	  if (qryMatchRefStart == -1) {
	    qryMatchRefStart = chainRows[0].tStart;
	  }
	  if (qryMatchRefEnd == -1) {
	    qryMatchRefEnd = chainRows[chainRows.size() -1 ].tEnd;
	  }
	  // Compare the score of the potential inversion with that of the chained seq.
	  chainScore = ScoreRows(tSeq, qSeq, chainRows, qryMatchRefStart, qryMatchRefEnd, 
				 lengthAligned, numMatches, numMismatches, scoreMat);
	  if (verbose > 1) {
	    std::cout << "query score on the chain, pos:" << qryMatchRefStart << " to: " << qryMatchRefEnd
		      << " to query : " << chainScore << " compare to inv: " 
		      << invScore << std::endl;
	    std::cout << "align stats: length aligned: " << lengthAligned << " matches: " 
		      << numMatches << " mismatches " << numMismatches 
		      << " compare to " << invLengthAligned << " "
		      << invNumMatches << " " << invNumMismatches << std::endl;
	    std::cout << "scores: " << invScore << " " << chainScore
		      << " " << chainScore  - (avgPenalty*0.5*(numMismatches))
		      << " allowing for " << (avgPenalty*0.5*(numMismatches)) << std::endl;
	  }
	  if (invScore >= (chainScore - (avgPenalty*0.5*(numMismatches )))
	      and lengthAligned *2 > invLengthAligned) {
	  //	  if (invScore >= chainScore) {
	    if (verbose) {
	      std::cout << "better score on the chain " << invScore << " " << chainScore
			<< " " << chainScore  - (avgPenalty*0.5*(numMismatches))
			<< " allowing for " << (avgPenalty*0.5*(numMismatches)) << std::endl;
	    }
	    if (explain) {
	      outFile << "EX: " << lavFileName << " " << refStart << " " << refEnd 
		      << " better alignment on query " << std::endl;
	    }
	    continue;
	  }
	}
	else {
	  if (verbose) {
	    std::cout << "didin't find any foward strand alignment with: "
		      << query.prevQuery << std::endl;
	  }
	}
	      
      
	LAVBlock *trueInversion = new LAVBlock;
	trueInversion->refBegin = refStart;
	trueInversion->refEnd   = refEnd;
	trueInversion->qryBegin = qryStart;
	trueInversion->qryEnd   = qryEnd;
	trueInversions.push_back(trueInversion);
      }
    }
  }
  ssize_t i;
  for (i = 0; i < trueInversions.size(); i++) {
    block = trueInversions[i];
    outFile << "OK: " << lavFileName << " " << block->refBegin << " " << block->refEnd << " " 
	    << block->qryBegin << " " << block->qryEnd << " " 
	    << orientation << std::endl;
    
  }
}

ssize_t GetInversionRows(DBQuery &query, std::string tableName, 
		     ssize_t refStart, ssize_t refEnd, 
		     ssize_t qryStart, ssize_t qryEnd,
		     std::vector<LAVRow> & inversion) {
  LAVRow lavRow;
  MYSQL_ROW row;
  MYSQL_RES *res;
  query << "SELECT * FROM " << tableName 
	<< " WHERE (qStart >= " << qryStart << " AND qEnd <= " << qryEnd << ") "
	<< " AND (tStart >= " << refStart << " AND tEnd <= " << refEnd << ") "
	<< " AND strand = 1 ";
  query.execute();
  res = mysql_store_result(query.conn);
  while (row = mysql_fetch_row(res)) {
    GetBlock(row, lavRow);
    inversion.push_back(lavRow);
    /*
      std::cout << "got inversion row: " << lavRow.tStart << " " << lavRow.tEnd << " " 
      << lavRow.qStart << " " << lavRow.qEnd << std::endl;
    */
  }
}

void GetOverlappingCoordinates(DBQuery &query, std::string &netTableName,
			       ssize_t refLength, ssize_t qryLength,
			       ssize_t &refStart, ssize_t &refEnd, ssize_t &qryStart, ssize_t &qryEnd) {

  query << "select max(score) from " << netTableName;
  ssize_t maxScore;
  GetOneInt(query, maxScore);
  query << "select * from " << netTableName << " where score= " << maxScore;
  query.execute(); 

  MYSQL_RES *res;
  MYSQL_ROW row;
  refStart = 0; refEnd = refLength; qryStart = 0; qryEnd = qryLength;
  
  if ((res = mysql_store_result(query.conn)) == NULL)
    return;

  row = mysql_fetch_row(res);
  if (row == NULL) 
    return;
  ssize_t netTStart, netTEnd, netQStart, netQEnd;
  netTStart = atoi(row[6]);
  netTEnd   = atoi(row[7]);
  netQStart = atoi(row[8]);
  netQEnd   = atoi(row[9]);
  
  // Now figure out how the two genomes map to eachother.
  if (netTStart < netQStart ) {
    refStart = 0;
    qryStart = netQStart - netTStart; //assume no gaps
  }
  else {
    qryStart = 0;
    refStart = netTStart - netQStart;
  }

  if ((refLength - netTEnd) < (qryLength - netQEnd)) {
    refEnd = refLength -1;
    qryEnd = netQEnd + (refLength - netTEnd);
  }
  else {
    qryEnd = qryLength - 1;
    refEnd = netTEnd + (qryLength - netQEnd);
  }
}
