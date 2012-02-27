/***************************************************************************
 * Title:          dbInvEnumerate.cpp 
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
  std::cout << "usage: dbinvcheck  [-d database -s seqNameTable " << std::endl
	    << "                    -w window]  lavFile." << std::endl;
}

void GetOverlappingCoordinates(DBQuery &query, std::string &netTableName,
			       ssize_t refLength, ssize_t qryLength,
			       ssize_t &refStart, ssize_t &refEnd, 
			       ssize_t &qryStart, ssize_t &qryEnd);


ssize_t InitEnv(int argc, char* argv[], 
	    std::string &dbName,
	    std::string &lavFile,
	    std::string &outFileName,
 	    std::string &seqTableName,
	    std::string &scoreMatFile,
	    ssize_t &window,
	    ssize_t &verbose) {
  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "d:s:w:m:v")) != EOF){
    switch(copt) {
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
  LAVRow lavBlock;
  dbName       = "july_alignments";
  seqTableName = "sequences";
  char *home = getenv("HOME");
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";

  ssize_t window = 100;
  ssize_t verbose = 0;
  InitEnv(argc, argv, 
	  dbName, lavFileName, outputFileName, seqTableName,
	  scoreMatName, 
	  window,
	  verbose);

  double scoreMat[5][5];
  double **scoreMatPtr;

  std::ofstream outFile;
  openck(outputFileName, outFile, std::ios::app);

  ReadScoreMatFile(scoreMatName, scoreMat);
  AssignScoreMatPtr(scoreMatPtr, scoreMat);
  
  ParseFileName(lavFileName, refSpec, qrySpec, chrom);
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

  std::string tFileName, qFileName;
  tFileName = refSpec + "." + chrom + ".fa";
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
  
  ssize_t foundInSelf;
  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);
  
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
  if (CreateTemporaryTable(query,
			   tableName, chainTableName, netTableName, tmpTableName,
                           0,
			   tempTables, 1, 0) == 0) {
    // better handled by exceptions, but that's a whole other bag of worms.
    // Don't output anything since that would imply that there are no 
    // inversions between these two sequences, instead, no comparison was possible.
    return 0;
  }
  
  
  if (verbose == 2) {
    std::cout << "creating " << tmpTableName << std::endl;
    std::cout << query.str() << std::endl;
  }

  //  query.execute();
  // Get the length of the query string
  /*
  query << "CREATE TEMPORARY TABLE " << refSelfTable << " AS SELECT " 
	<< refSelfTable << ".id, " 
	<< refSelfTable << ".tStart, " 
	<< refSelfTable << ".tEnd, "
	<< refSelfTable << ".qStart, "
	<< refSelfTable << ".qEnd, "
	<< refSelfTable << ".strand, "
	<< refSelfChain << ".chainId " 
	<< " FROM " << refSelfTable << ", " << refSelfChain 
	<< " WHERE " << refSelfTable << ".id = " << refSelfChain << ".lavid AND "
	<< refSelfTable << ".strand = 1 ";
  if (verbose) {
    std::cout << "creating table with " << query.str() << std::endl;
  }
  query.execute();

  query << "CREATE TEMPORARY TABLE " << qrySelfTable << " AS SELECT " 
	<< qrySelfTable << ".id, " 
	<< qrySelfTable << ".tStart, " 
	<< qrySelfTable << ".tEnd, "
	<< qrySelfTable << ".qStart, "
	<< qrySelfTable << ".qEnd, "
	<< qrySelfTable << ".strand," 
	<< qrySelfChain << ".chainId " 
	<< " FROM " << qrySelfTable << ", " << qrySelfChain 
	<< " WHERE " << qrySelfTable << ".id = " << qrySelfChain << ".lavid AND "
	<< qrySelfTable << ".strand = 1 ";
  if (verbose) {
    std::cout << "creating qry self chain: " << query.str() << std::endl;
  }
  query.execute();
  */
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
		 ( (block2->qryBegin - qryEnd) < 500) and
		 ( (block2->refBegin - refEnd) < 500) ) {
	      a = a2; // skip considering this inversion.
	      if (verbose ){
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
	  for (a2 = a+1; a2 < alignedContig->alignments.size(); a2++) {
	    block2 = alignedContig->alignments[a2];
	    if (block2->refEnd == netTEnd) {
	      a = a2+1;
	      break;
	    }
	  }
	}

	// The coordinates of the ref self-alignment 
	// are always the same (regardless of the orientation of the
	// qry). 
      
	if (verbose) {
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
	    std::cout << "skipping because coordinates are outside of orthologous region " << std::endl;
	    std::cout << refStart << " " << refEnd << " " << qryForwardStart 
		      << " " << qryForwardEnd << std::endl;
	    std::cout << refOvpStart << " " << refOvpEnd << " " << qryOvpStart 
		      << " " << qryOvpEnd << std::endl;
	  }
	  continue;
	}
	    
	std::vector<LAVRow> inversionRows;
	GetInversionRows(query, tableName, refStart, refEnd, qryStart, qryEnd, inversionRows);
	invScore = ScoreRows(tSeq, qSeq, inversionRows, -1, -1, scoreMatPtr);

	// Get all of the rows corresponding to the inversion.

	// Check 1.  
	// Look to see if there is an alignment in the reference
	// sequence that covers the same area as the inversion. If
	// this is the case, the inversion is just a palindrome. 

	// The tmprefselftable only has the reverse strand alignments

	// It's possible that many chained alignments are present at the 
	// position of the inversion. 


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
	      palindromeScore = ScoreRows(tSeq, tSeq, palindromeRows, -1, -1, scoreMatPtr);
	      if (verbose == 2) {
		std::cout << "got ref palindrome score: " << palindromeScore << " inversion score: " 
			  << invScore<< std::endl;
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
	  continue;
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
	      palindromeScore = ScoreRows(qSeq, qSeq, palindromeRows, -1, -1, scoreMatPtr);
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
	  continue;
	}

	mysql_free_result(res);
      
      
	// Now look to see if any alignments happen in the same row (appear twice 
	// in the query sequence).
	// Assume that the query alignment in questin is in the reverse strand (although
	// this isn't the case for sequenes that are already in oppposite direction).
      
	// Look to see if this inversion is on the diagonal.

	// Get the beginning of the alignment that is closest to the 
	// middle of the inversion.
	LAVRow prevBlock, nextBlock;
	ssize_t maxTStart;

	// Now check to see if there is a sequence corresponding to a chain that 
	// covers the position of the inversion.

	query << "SELECT " << tmpTableName << ".id, tStart, tEnd, qStart, qEnd, strand, chainId" 
	      << " from " << tmpTableName 
	      << " WHERE ((( tStart >= " << refStart 
	      << " AND tStart <= " << refEnd << ")"
	      << "        OR ( tEnd >= " << refStart 
	      << " AND tEnd <= " << refEnd << ") "
	      << "        OR ( tStart <= " << refStart 
	      << " AND tEnd >= "  << refEnd << "))  "
	      << "     OR (( qStart >= " << qryForwardStart 
	      << " AND qStart <= "  << qryForwardEnd << ")"
	      << "      OR ( qEnd >= " << qryForwardStart 
	      << " AND qEnd <= " << qryForwardEnd << ")"
	      << "      OR ( qStart <= " << qryForwardStart 
	      << " AND qEnd >= " << qryForwardEnd << ")))"
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
	  if (verbose == 2) {
	    std::cout << "got row to compare: " << lavRow.tStart << " " << lavRow.tEnd << " "
		      << lavRow.qStart << " " << lavRow.qEnd << std::endl;
	  }
	}
	double  chainScore;
	if (chainRows.size() > 0) {
	  // Compare the score of the potential inversion with that of the chained seq.
	  chainScore = ScoreRows(tSeq, qSeq, chainRows, refStart, refEnd, scoreMatPtr);
	  if (invScore >= (chainScore - 1000)) {
	    if (verbose) {
	      std::cout << "better score on the chain " << invScore << " " << chainScore << std::endl;
	    }
	    continue;
	  }
	}
      
	LAVBlock *trueInversion = new LAVBlock;
	trueInversion->refBegin = refStart;
	trueInversion->refEnd   = refEnd;
	trueInversion->qryBegin = qryStart;
	trueInversion->qryEnd   = qryEnd;
	//	trueInversion->enumeration = enumeration;
	trueInversions.push_back(trueInversion);
	//	++enumeration;
      }
    }
  }
  ssize_t i;
  outFile << lavFileName << std::endl;
  for (i = 0; i < trueInversions.size(); i++) {
    block = trueInversions[i];
    outFile << block->refBegin << " " << block->refEnd << " " 
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

