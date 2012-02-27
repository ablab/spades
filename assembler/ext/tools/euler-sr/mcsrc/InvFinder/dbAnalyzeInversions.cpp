/***************************************************************************
 * Title:          dbAnalyzeInversions.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <stdlib.h>

// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <algorithm>

#include <unistd.h>

#include <ctype.h>

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
#include "blocks/BlockDB.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "TupleLib.h"
#include "alignutils.h"
#include "SeqUtils.h"
#include "emboss/EmbossAlignment.h"
#include "emboss/EmbossAlign.h"

ssize_t dbDetermineOrientation(DBQuery &query, std::string &tableName) {
  ssize_t plus;
  ssize_t minus;
  query << "SELECT COUNT(chainId) FROM "<< tableName << " WHERE strand=0";
  GetOneInt(query, plus);
  query << "SELECT COUNT(chainId) FROM "<< tableName << " WHERE strand=1";
  GetOneInt(query, minus);

  if (plus > minus) return 0;
  else return 1;
}



class DBDescription {
public:
  std::string dbName;
  std::map<std::string, std::string> specTable;
};

typedef std::map<std::string, DBDescription*> DBDescriptionMap;

void ReadInversionFile(std::string fileName, 
		       std::vector<std::string> &refNames,
		       std::vector<std::string> &qryNames,
		       std::vector<std::string> &seqNames,
		       std::vector<ssize_t> &numInversions,
		       std::vector<ssize_t> &refStartPositions,
		       std::vector<ssize_t> &refEndPositions,
		       std::vector<ssize_t> &qryStartPositions,		       
		       std::vector<ssize_t> &qryEndPositions);
void PrintUsage() {
  std::cout << "dbinvcheck: database based inversion analysis." << std::endl;
  std::cout << "use this to check the conservation surrounding inversions. " << std::endl;
  std::cout << "usage: dbai  [-d database -s seqNameTable -B borerlen -w window  " 
	    <<  "-a firstout  -b secondout (fragment files) " << std::endl;
  std::cout << " -B borderlen, length of seq at borders of alignment to align " << std::endl;
  std::cout << " -g gap open penalty (400) " << std::endl;
  std::cout << " -e gap extend penalty (30) " << std::endl;
  std::cout << " -m scoreMatrix file to use with score matrix " << std::endl;
  std::cout << " -w window of scores to compute average of. " << std::endl;
}

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &dbName,
	    std::string &inversionFileName,
	    std::string &seqTableName,
	    std::string &firstOut,
	    std::string &secondOut,
	    std::string &scoreMatFile,
	    double &gapOpen,
	    double &gapExtend,
	    ssize_t &borderLen,
	    ssize_t &window) {
  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "o:d:s:m:w:a:b:g:e:B:")) != EOF){
    switch(copt) {
    case 'B':
      borderLen = atoi(optarg);
      continue;
    case 'g':
      gapOpen = atof(optarg);
      continue;
    case 'e':
      gapExtend = atof(optarg);
      continue;
    case 'a': 
      firstOut = optarg;
      continue;
    case 'b':
      secondOut = optarg;
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
    case 'm':
      scoreMatFile = optarg;
      continue;
    default:
      PrintUsage();
      exit(0);
      continue;
    }
  }
  i = optind;
  if (i < argc) { 
    inversionFileName = argv[i++];
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



int main(int argc, char* argv[]) {
  std::string outputFileName;
  std::string seqTableName;
  std::string dbName;
  std::string invFileName;
  std::string invScoreFileName;
  std::string refSpec, qrySpec, chrom;
  
  dbName       = "lav_alignments";
  seqTableName = "sequences";

  FloatMatrix scoreMat;
  ssize_t window = 50;
  ssize_t borderLen = 1000;
  std::vector<std::string> refNames, qryNames, seqNames;
  std::vector<ssize_t> numInversions, refStartPos, refEndPos, qryStartPos, qryEndPos;

  std::string outA, outB;
  outA = "first.fasta";
  outB = "second.fasta";
  
  char *home;
  std::string scoreMatFileName = "";
  double gapOpen   = 400;
  double gapExtend = 30;

  home = getenv("HOME");
  if (home != NULL) 
    scoreMatFileName = std::string(home) 
      + "/projects/src/align/data/scoremat.txt";
  
  std::cout << "default score matrix: " << scoreMatFileName << std::endl;
  InitEnv(argc, argv, 
	  dbName, invFileName, seqTableName, 
	  outA, outB,
	  scoreMatFileName,
	  gapOpen, gapExtend,
	  borderLen, window);

  invScoreFileName = invFileName + ".scores";
  std::ofstream invScoreFile;
  openck(invScoreFileName, invScoreFile);

  if (scoreMatFileName != "" ) {
    ReadScoreMatFile(scoreMatFileName, scoreMat);
  }

  ReadInversionFile(invFileName, refNames, qryNames, seqNames,
		    numInversions, refStartPos, refEndPos, qryStartPos, qryEndPos);

  ssize_t s, p,i;
  /* 
     Sanity check.  Print the inversion file if things are strange .
     i = 0;
     for (s = 0; s < refNames.size(); s++) {
     std::cout << refNames[s] << " " << qryNames[s] << " " << seqNames[s] << " " << numInversions[s];
     for (p = 0; p < numInversions[s] ; p++) {
     std::cout << refStartPos[i] << " " << refEndPos[i] << " " << qryStartPos[i] << " " << qryEndPos[i] << std::endl;
     i++;
     }
     }
  */

  
  // Locate the chains that make up each inversion.
  
  ssize_t foundInSelf;
  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);
  std::string tableName, chainTableName, tmpTableName;
  // Create the query and associate it with the databae
  DBQuery query(dbConn);

  MYSQL_RES *res;
  MYSQL_ROW row;
  ssize_t numRows;
  i = 0;
  ssize_t t;
  std::vector<ssize_t> chainIds, refStart, refEnd, qryStart, qryEnd;
  std::vector<LAVRow> lavRows;
  LAVRow lavRow;

  std::string refSeqName, qrySeqName;
  DNASequence refSeq, qrySeq, refInv, qryInv;

  refInv._ascii = 1;
  qryInv._ascii = 1;
  std::ofstream file1, file2;
  openck(outA, file1);
  openck(outB, file2);
  EmbossAlignment alignment;

  for (s = 0; s < refNames.size(); s++) {
    std::cout << "------------------------------------------------------------" << std::endl;
    std::cout << refNames[s] << " " << qryNames[s] 
	      << " " << seqNames[s] << " " << numInversions[s] << std::endl;
    tableName = refNames[s] + "_" + qryNames[s] + "_" + seqNames[s];
    chainTableName = tableName + "_chain";
    tmpTableName = tableName + "_tmp";
    // Load the reference and query sequences    

    // find out what lav alignments are in the chain table
    query << "CREATE TEMPORARY TABLE " << tmpTableName << " AS SELECT " 
	  << tableName << ".id, tStart, tEnd, qStart, qEnd, strand, " 
	  << tableName << ".chainId FROM " << tableName << ", " 
	  << chainTableName << " WHERE " << tableName << ".chainId = " 
	  << chainTableName << ".lavid";
    //  std::cout << "creating " << tmpTableName << std::endl;
    query.execute();

    refSeqName = refNames[s] + "." + seqNames[s] + ".fa"; 
    qrySeqName = qryNames[s] + "." + seqNames[s] + ".fa"; 

    SeqReader::GetSeq(refSeqName, refSeq, SeqReader::noConvert); 
    SeqReader::GetSeq(qrySeqName, qrySeq, SeqReader::noConvert); 
    
    for (p = 0; p < numInversions[s] ; p++) {
      
      // Things to work with here: 
      // ref{Start,End}Pos, qry{Start,End}Pos
      /*
	refInv.seq = &refSeq.seq[refStartPos[i]];
	refInv.length = refEndPos[i] - refStartPos[i] + 1;

	qrytInv.seq = &qrySeq.seq[qrySeq.length - qryStartPos[i]];
	qryInv.length = qryEndPos[i] - qryStartPos[i] + 1;
      */
      query << "SELECT * FROM " << tableName << " WHERE ((qStart >= " << qryStartPos[i] 
	    << " AND qEnd <= " << qryEndPos[i] << ") OR ( qStart <= " << qryEndPos[i] 
	    << " AND qEnd >= " << qryEndPos[i] << " ) OR ( qStart <= " 
	    << qryEndPos[i] << " AND qEnd > " << qryEndPos[i] << ")) AND strand=1";

      query.execute();
      res = mysql_store_result(query.conn);
      numRows = mysql_num_rows(res);
      ssize_t r;
      for (r = 0; r < numRows; r++) {
	row = mysql_fetch_row(res);
	chainIds.push_back(atoi(row[0]));
	GetBlock(row, lavRow);
	lavRows.push_back(lavRow);
      }

      // Find the chains that end just before and begin just after this repeat.

      ssize_t prevChainId, nextChainId;
      ssize_t tStart, tEnd;
      if (chainIds.size() > 0) {
	prevChainId = chainIds[0]-1;
	nextChainId = chainIds[chainIds.size()-1]+1;

	LAVRow prevChainRow, nextChainRow;
	ssize_t id;

	// Get the blocks ending just before and starting just after the inversion.
	query << "SELECT max(tEnd) FROM " << tmpTableName 
	      << " WHERE tEnd < " << lavRows[0].tStart << " AND strand=0";
	GetOneInt(query, tEnd);
	query << "SELECT * FROM " << tableName << " WHERE tEnd = " << tEnd;
	res = query.executeNF();
	row = mysql_fetch_row(res);
	GetBlock(row, prevChainRow);

	query << "SELECT MIN(tStart) FROM " << tableName 
	      << " WHERE tStart >  " << lavRows[lavRows.size()-1].tEnd << " AND strand=0";
	GetOneInt(query, tStart);
	query << "SELECT * FROM  " << tableName << " WHERE tStart = " << tStart << " AND strand=0";
	res = query.executeNF();
	row = mysql_fetch_row(res);
	GetBlock(row, nextChainRow);

	// Now check to see how long the gaps are between alignments.

	ssize_t prevRefGap, nextRefGap, prevQryGap,  nextQryGap;
	ssize_t qrySeqLength;

	prevRefGap = lavRows[0].tStart - prevChainRow.tEnd + 1;
	nextRefGap = nextChainRow.tStart - lavRows[lavRows.size()-1].tEnd + 1;
	
	
	query << "SELECT length from " << seqTableName << " WHERE name=" 
	      << dbstr(qryNames[s]) << " AND sequence=" << dbstr(seqNames[0]);
	ssize_t queryLength;
	GetOneInt(query, queryLength);
	ssize_t orientation;
	orientation = dbDetermineOrientation(query, tableName);

	// if (orientation == 0) {
	// query strand is in the same orientation as the reference strand, so the
	// inversions are in the opposite direction.  They then need to be mapped to
	// the other direction to compare their start positiion to the previous block;
	ssize_t invQryStart, invQryEnd; // The coordinates of the inversion on the forward strand.
	// todo: make this work for sequences that are in the reverse direction.
	invQryStart = queryLength - qryEndPos[i];
	invQryEnd   = queryLength - qryStartPos[i];
	
	prevQryGap = invQryStart - prevChainRow.qEnd + 1;
	nextQryGap = nextChainRow.qStart - invQryEnd + 1; 

	// Query strand is in the revsrse orientation of the reference strand
	// the inversions are in the same dire ction.
	// But nothing is changed about the lengths of the gaps between alignments.
	
	// Do all sorts of alignments
	std::cout << "Checking out sequence:" 
		  << " " << refStartPos[i] << " " << refEndPos[i] 
		  << " " << qryStartPos[i] << " " << qryEndPos[i] 
		  << " " << prevRefGap << " " << prevQryGap << std::endl;
	
	double alignScore;
	double pctIdentity, repPctIdentity, uniqPctIdentity;
	double refRepeatContent, qryRepeatContent;
	ssize_t *alignLocations = NULL;
	ssize_t lastRow;
	lastRow = lavRows.size()-1;
	DNASequence refFrag, qryFrag;
	refFrag._ascii = 1;
	qryFrag._ascii = 1;
#ifdef RUN_ALIGNMENTS
	if (prevRefGap > 0 && prevRefGap < 10000 && 
	    prevQryGap > 0 && prevQryGap < 10000) {

	  ssize_t refStart, qryStart;
	  // Output the aligned block before the gap before the inversion
	  if (prevChainRow.tEnd - prevChainRow.tStart+1 < borderLen) {
	    refStart = std::max(0, prevChainRow.tEnd - borderLen);
	    refFrag.seq = &refSeq.seq[refStart];
	    refFrag.length = prevChainRow.tEnd - refStart+ 1;
	  }
	  else {
	    refStart = prevChainRow.tStart;
	    refFrag.seq = &refSeq.seq[prevChainRow.tStart];
	    refFrag.length = prevChainRow.tEnd - prevChainRow.tStart + 1;
	  }
	  *refFrag.titlestream << tableName << "_A_" << refStart << "_" << prevChainRow.tEnd;
	  refFrag.PrintSeq(file1);
	  file1 << std::endl;
	
	  // Output the query sequence of the aligned block before the inversion
	  if (prevChainRow.qEnd - prevChainRow.qStart+1 < borderLen) {
	    qryStart = std::max(0, prevChainRow.qEnd - borderLen);
	    qryFrag.seq = &qrySeq.seq[qryStart];
	    qryFrag.length = prevChainRow.qEnd - qryStart;
	  }
	  else {
	    qryFrag.seq = &qrySeq.seq[prevChainRow.qStart];
	    qryFrag.length = prevChainRow.qEnd - prevChainRow.qStart + 1;
	  }
	  *qryFrag.titlestream << tableName << "_A_" << qryStart << "_" << prevChainRow.qStart;
	  qryFrag.PrintSeq(file2);
	  file2 << std::endl;
	  alignLocations = NULL;
	  if (refFrag.length * qryFrag.length < 1000000) {
	    alignScore =  AffineAlign(refFrag, qryFrag, 
				      -100, 100, 200,
				      gapOpen, gapExtend,
				      alignLocations, scoreMat);
	
	    pctIdentity = CalcPercentIdentity(refFrag, qryFrag, alignLocations);
	    delete []alignLocations;
	    alignLocations = NULL;

	    std::cout << " block/forward: "
		      << refStart << "-" << prevChainRow.tEnd  
		      << "(" << prevChainRow.tEnd - refStart + 1 << ")"
		      << ", " << qryStart << "-" << prevChainRow.qEnd  
		      << "(" << prevChainRow.qEnd - qryStart << ")" 
		      << " score: " << alignScore 
		      << " pctidentity: " << pctIdentity
		      << std::endl;
	  }
	  else {
	    std::cout << " block/forward: "
		      << refStart << "-" << prevChainRow.tEnd  
		      << "(" << prevChainRow.tEnd - refStart + 1 << ")"
		      << ", " << qryStart << "-" << prevChainRow.qEnd  
		      << "(" << prevChainRow.qEnd - qryStart << ")" 
		      << std::endl;
	  }

	  // Output the gapped sequence before the inversion
	  refFrag.seq = &refSeq.seq[prevChainRow.tEnd];
	  refFrag.length = prevRefGap;
	  *refFrag.titlestream << tableName << "_B_" << refStartPos[i];
	  refFrag.PrintSeq(file1);
	  file1 << std::endl;

	  // Output the query sequence of the unaligned gap before the inversion
	  qryFrag.seq = &qrySeq.seq[prevChainRow.qEnd];
	  qryFrag.length = prevQryGap;
	  *qryFrag.titlestream << tableName << "_B_" << qryStartPos[i];
	  qryFrag.PrintSeq(file2);
	  file2 << std::endl;
	  if (refFrag.length * qryFrag.length < 1000000) {
	    alignLocations = NULL;
	    alignScore =  AffineAlign(refFrag, qryFrag,  
				      -100, 100, 200,
				      gapOpen, gapExtend,
				      alignLocations, scoreMat);
	  
	    pctIdentity = CalcPercentIdentity(refFrag, qryFrag, alignLocations);
	    delete []alignLocations;
	    alignLocations = NULL;
	  
	    std::cout << " forward/inverted: " << prevChainRow.tEnd << "-" << refStartPos[i] 
		      << "(" << prevRefGap << ")"
		      << ", " << prevChainRow.qEnd << "-" << qryStartPos[i]
		      << "(" << prevQryGap << ")"
		      << " score: " << alignScore 
		      << " pctidentity: " << pctIdentity << std::endl;
	  }
	  else {
	    std::cout << " forward/inverted: " << prevChainRow.tEnd << "-" << refStartPos[i] 
		      << "(" << prevRefGap << ")"
		      << ", " << prevChainRow.qEnd << "-" << qryStartPos[i]
		      << "(" << prevQryGap << ")"
		      << std::endl;
	  }
	}
      
	// Output the ref sequence corresponding to the inversion
	refFrag.seq = &refSeq.seq[refStartPos[i]];
	refFrag.length = refEndPos[i] - refStartPos[i];
	*refFrag.titlestream << tableName << "_C_" << lavRows[0].tStart 
			    << "_" << lavRows[lastRow].tEnd ;
	refFrag.PrintSeq(file1);
	file1 << std::endl;

	// Output the qry sequence corresponding to the inversion
	qryFrag.seq = &qrySeq.seq[qrySeq.length - qryEndPos[i]];
	qryFrag.length = qryEndPos[i] - qryStartPos[i] + 1;
	MakeRC(qryFrag, qryInv);
	*qryInv.titlestream << tableName << "_C_" << qryStartPos[i] 
			   << "_" << qryEndPos[i] ;
	qryInv.PrintSeq(file2);
	file2 << std::endl;
	if (refFrag.length * qryInv.length < 10000000) {
	  alignScore =  AffineAlign(refFrag, qryInv,  
				    -100, 100, 200,
				    gapOpen, gapExtend,
				    alignLocations, scoreMat);
	  refRepeatContent = CountRepeatMasked(refFrag);
	  qryRepeatContent = CountRepeatMasked(refInv);
	  pctIdentity = CalcPercentIdentity(refFrag, qryInv, alignLocations);
	  delete []alignLocations;
	  alignLocations = NULL;
	  std::cout << " inversion: " << lavRows[0].tStart << "-" << lavRows[lastRow].tEnd
		    << "(" << lavRows[lastRow].tEnd - lavRows[0].tStart << ")"
		    << ", " << lavRows[0].qStart << " " << lavRows[lastRow].qEnd  
		    << "(" << lavRows[lastRow].qEnd - lavRows[0].qStart << ")"
		    << " score: " << alignScore 
		    << " pctidentity: " << pctIdentity
		    << std::endl;
	}
	else {
	  std::cout << " inversion: " << lavRows[0].tStart << "-" << lavRows[lastRow].tEnd
		    << "(" << lavRows[lastRow].tEnd - lavRows[0].tStart << ")"
		    << ", " << lavRows[0].qStart << " " << lavRows[lastRow].qEnd  
		    << "(" << lavRows[lastRow].qEnd - lavRows[0].qStart << ")"
		    << std::endl;
	}
      
      if (nextRefGap > 0 && nextRefGap < 10000 &&
	  nextQryGap > 0 && nextQryGap < 10000) {

	// Output the sequence that is a gap in the alignment after the 
	// inversion.
	refFrag.seq = &refSeq.seq[lavRows[lastRow].tEnd];
	refFrag.length = nextRefGap;
	*refFrag.titlestream << tableName << "_D_"<< lavRows[lastRow].tEnd;
	refFrag.PrintSeq(file1);
	file1 << std::endl; 
	  
	// Output the query sequence of the unaligned region after 
	// the inversion
	qryFrag.seq = &qrySeq.seq[lavRows[lastRow].qEnd];
	qryFrag.length = nextQryGap;
	*qryFrag.titlestream << tableName << "_D_" << lavRows[lastRow].qEnd;
	qryFrag.PrintSeq(file2);
	file2 << std::endl;
	  
	if (refFrag.length * qryFrag.length < 10000000) {
	  alignScore =  AffineAlign(refFrag, qryFrag,  
				    -100, 100, 200,
				    gapOpen, gapExtend,
				    alignLocations, scoreMat);
	  pctIdentity = CalcPercentIdentity(refFrag, qryFrag, alignLocations);
	  delete []alignLocations;
	  alignLocations = NULL;
	  std::cout << " inverted/next: " << lavRows[lastRow].tEnd << "-" << nextChainRow.tStart
		    << "(" << nextChainRow.tStart - lavRows[lastRow].tEnd << ")"
		    << ", " << lavRows[lastRow].qEnd + 1 << "-" << nextChainRow.qStart
		    << "(" << nextChainRow.qStart - lavRows[lastRow].qEnd + 1 << ")"
		    << " score: " << alignScore 
		    << " pctidentity: " << pctIdentity
		    << std::endl; 
	}
	else {
	  std::cout << " inverted/next: " << lavRows[lastRow].tEnd << "-" << nextChainRow.tStart
		    << "(" << nextChainRow.tStart -  lavRows[lastRow].tEnd << ")"
		    << ", " << lavRows[lastRow].qEnd << "-" << nextChainRow.qStart  
		    << "(" << nextChainRow.qStart - lavRows[lastRow].qEnd + 1 << ")"
		    << std::endl;
	}
	
	ssize_t qryEnd, refEnd;
	// Output the block of the alignment after the inversion
	if (nextChainRow.tStart + borderLen > nextChainRow.tEnd)
	  refEnd = std::min(nextChainRow.tStart+borderLen -1, refSeq.length-1);
	else
	  refEnd = nextChainRow.tEnd;
      
	if (nextChainRow.qStart + borderLen > nextChainRow.tEnd)
	  qryEnd = std::min(nextChainRow.qStart + borderLen - 1, qrySeq.length-1);
	else
	  qryEnd = nextChainRow.qEnd;

	refFrag.seq = &refSeq.seq[nextChainRow.tStart];
	refFrag.length = refEnd - nextChainRow.tStart + 1;
	*refFrag.titlestream << tableName << "_E_" << nextChainRow.tStart;
	refFrag.PrintSeq(file1);
	file1 << std::endl;

	// Output the block of the qry seq the alignment after the inversion
	qryFrag.seq = &qrySeq.seq[nextChainRow.qStart];
	qryFrag.length = qryEnd - nextChainRow.qStart + 1;
	*qryFrag.titlestream << tableName << "_E_" << nextChainRow.qStart;
	qryFrag.PrintSeq(file2);
	file2 << std::endl;
	if (refFrag.length * qryFrag.length < 10000000) {
	  alignScore =  AffineAlign(refFrag, qryFrag,  
				    -100, 100, 200,
				    gapOpen, gapExtend,
				    alignLocations, scoreMat);
	  pctIdentity = CalcPercentIdentity(refFrag, qryFrag, alignLocations);
	  delete []alignLocations;
	  alignLocations = NULL;
	  std::cout << " last block: " << nextChainRow.tStart << " - " << refEnd
		    << "(" << refEnd - nextChainRow.tStart << ")"
		    << ", " << nextChainRow.qStart + 1 << " - " << qryEnd
		    << "(" << qryEnd - nextChainRow.qStart + 1 << ")"
		    << " score: " << alignScore 
		    << " pctidentity: " << pctIdentity
		    << std::endl; 
	}
	else {
	  std::cout << " last block: " << nextChainRow.tStart << " - " << refEnd
		    << "(" << refEnd - nextChainRow.tStart << ")"
		    << ", " << nextChainRow.qStart + 1 << " - " << qryEnd
		    << "(" << qryEnd - nextChainRow.qStart + 1 << ")"
		    << std::endl;
	}
      }  // End checking to see if the alignment size is too large.
      #endif
      // Try to do the full alignment, except with the inversion in the proper orientation.
      ssize_t refStart, refEnd;
      ssize_t qryStart, qryEnd;

      /*
	if (prevChainRow.tStart +200 > prevChainRow.tEnd)
	refStart = prevChainRow.tStart;
      else
      */
      refStart = prevChainRow.tEnd - borderLen;

      /*
	if (prevChainRow.qStart + borderLen > prevChainRow.qEnd)
	qryStart = prevChainRow.qStart;
      else
      */
	qryStart = prevChainRow.qEnd - borderLen;

      /*
	if (nextChainRow.tStart + borderLen > nextChainRow.tEnd)
	refEnd = nextChainRow.tEnd;
      else
      */
      refEnd = nextChainRow.tStart + borderLen;

      /*      if (nextChainRow.qStart + borderLen < nextChainRow.qEnd)
	qryEnd= nextChainRow.qEnd;
      else
      */
      qryEnd = nextChainRow.qStart + borderLen;

      refFrag.seq = &refSeq.seq[refStart];
      refFrag.length = refEnd-  refStart + 1;

      DNASequence qryInvRev;
      qryFrag.seq = &qrySeq.seq[qryStart];
      qryFrag.length = qryEnd - qryStart + 1;

      // Now revert the inversion.
      qryInv.seq = &qrySeq.seq[invQryStart];
      qryInv.length = invQryEnd - invQryStart + 1;

      MakeRC(qryInv, qryInvRev);
      ssize_t l;
      std::ofstream tmpout;
      tmpout.open("queryBefore.fasta");
      qryFrag.PrintSeq(tmpout);
      tmpout << std::endl;
      tmpout.close();
      tmpout.clear();
      tmpout.open("inv.fasta");
      qryInv.PrintSeq(tmpout);
      tmpout<< std::endl;
      tmpout.close();

      for (l = 0; l < qryInv.length; l++)
	qryFrag.seq[l+invQryStart - qryStart] = qryInvRev.seq[l];

      if (alignment.locations != NULL)
	delete[] alignment.locations;

      std::string alignCommand;
      *refFrag.titlestream << "reference";
      *qryFrag.titlestream << "query";
      alignCommand = "~/projects/sw/bin/stretcher -aformat srspair ";
      ssize_t markers[4];
      markers[0] = prevChainRow.tEnd - refStart;
      markers[1] = refStartPos[i] - refStart;
      markers[2] = refEndPos[i] - refStart;
      markers[3] = nextChainRow.tStart - refStart;
      std::cout << markers[0] << " " << markers[1] << " " << markers[2] << " " << markers[3] << std::endl;
      if (!EmbossAlign(alignCommand, refFrag, qryFrag, alignment) ) {
	std::cout <<" error aligning sequence " << std::endl;
      }
      else {
	double *averageScores;
	ssize_t scoreLength;
	alignment.CalculateAverageScore(window, averageScores, scoreLength, 
					scoreMat, 300, 10, refFrag, qryFrag,
					markers, 4);

	ssize_t z;
	invScoreFile << scoreLength 
		     << " " << markers[0]
		     << " " << markers[1]
		     << " " << markers[2]
		     << " " << markers[3];
	
	for (z = 0; z < scoreLength; z++) 
	  invScoreFile << " " << averageScores[z];
	invScoreFile << std::endl;

	delete averageScores;
      }
      std::cout << "got alignemnt with score: " << alignment.alignScore << std::endl;
      // Align the sequences corresponding to the gaps.
      }
      lavRows.clear();
      chainIds.clear();
      std::cout << std::endl;
      i++;
    }
  }
  file1.close();
  file2.close();
  invScoreFile.close();
  std::cout << std::endl;
}

void ReadInversionFile(std::string fileName, 
		       std::vector<std::string> &refNames,
		       std::vector<std::string> &qryNames,
		       std::vector<std::string> &seqNames,
		       std::vector<ssize_t> &numInversions,
		       std::vector<ssize_t> &refStartPositions,
		       std::vector<ssize_t> &refEndPositions,
		       std::vector<ssize_t> &qryStartPositions,		       
		       std::vector<ssize_t> &qryEndPositions) {
  std::ifstream in;
  openck(fileName, in);
  std::string title;
  std::string refSpec, qrySpec, chrom;
  ssize_t ni; // num inversions
  ssize_t refStart, refEnd, qryStart, qryEnd;
  ssize_t strand;
  ni = -1;
  while (in) {
    if (isalpha(in.peek())) {
      // reading a title line.
      if ( !(in >> title)) {
	break;
      }
      if (ni != -1)
	numInversions.push_back(ni);
      ParseFileName(title, refSpec, qrySpec, chrom);      
      refNames.push_back(refSpec);
      qryNames.push_back(qrySpec);
      seqNames.push_back(chrom);
      ni = 0;
    }
    else {
      if (ni == -1) {
	std::cout << "got non-title information before the title, exiting "
		  << std::endl;
	exit(0);
      }
      if ( (in >> refStart >> refEnd >> qryStart >> qryEnd >> strand)) {
	ni++;
	refStartPositions.push_back(refStart);
	refEndPositions.push_back(refEnd);
	qryStartPositions.push_back(qryStart);
	qryEndPositions.push_back(qryEnd);
	in.get(); // discard the newline.
      }
    }
  }
  // add an entry for the last num inversions.
  if (ni != -1)
    numInversions.push_back(ni);
  in.close();
}
