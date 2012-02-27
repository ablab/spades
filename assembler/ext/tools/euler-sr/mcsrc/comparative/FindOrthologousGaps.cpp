/***************************************************************************
 * Title:          FindOrthologousGaps.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <unistd.h>

#include "DNASequence.h"
#include "SeqReader.h"
#include "blocks/BlockDB.h"
#include "mysql/mysql.h"

void GetTableName(std::string speciesA, std::string speciesB, 
		  std::string seqName,
		  std::string &tableName) {
  tableName = speciesA + "_" + speciesB + "_" + seqName;
}

void InitEnv(int argc, char* argv[], 
	     std::string &dbName,
	     std::string &refSpecies,
	     std::string &qrySpecies,
	     std::string &sequence,
	     std::vector<std::string> &species,
	     ssize_t &minGapLength,
	     ssize_t &maxGapLength,
	     ssize_t &minQryGapLength,
	     ssize_t &maxQryGapLength,
	     std::string &qryTableName);

void PrintUsage();

int main(int argc, char* argv[]) {
  // Steps to finding orthologous gaps.
  // 1. Take an alignment of sequences A and B, find a gap.
  // 2. Look at the orthologous positions of A in a sequence C.
  // 3. Find where C is aligned to B.  
  // 4. If the C has a gap relative to B that is of the same 
  //    or similar length as the gap to A, consider that an orthologous gap.

  ssize_t i; 
  i = 1;
  std::string dbName;
  std::string refSpecies;
  std::string qrySpecies;
  std::string seqName;
  ssize_t minGapLength, maxGapLength;
  std::vector<std::string> species;
  std::vector<DBQuery*> speciesQuery;
  std::string seqTableName;
  ssize_t minQryGapLength;
  ssize_t maxQryGapLength;

  dbName = "lav_alignments";
  seqTableName = "sequences";
  minGapLength = 100;
  maxGapLength = 10000;
  InitEnv(argc, argv, dbName, refSpecies, qrySpecies, seqName, 
	  species, 
	  minGapLength, maxGapLength, minQryGapLength, maxQryGapLength,
	  seqTableName);

  // Determine the names of the tables;
  std::string refVsQryTable, tableName;
  GetTableName(refSpecies, qrySpecies, seqName, refVsQryTable);

  MYSQL *conn;
  ConnectToDB(conn, dbName);
  
  std::cout << "fetching all alignments from: " << refVsQryTable << std::endl;
  DBQuery query(conn);
  query << "select * from " << refVsQryTable << " ORDER BY tStart " << std::endl;
  query.execute(); 
  
  MYSQL_RES *result;
  if ((result = mysql_store_result(conn)) == NULL) {
    std::cout << "error fetching result: " << mysql_error(conn) << " " << mysql_errno(conn) << std::endl;
  }
  ssize_t numRows;
  numRows = mysql_num_rows(result);
  std::cout << "done fetching " << numRows << " alignments from: " << refVsQryTable << std::endl;

  // Look for gaps in the reference vs query.

  // Step 1. Find out how many blocks to look through
  ssize_t numBlocks;
  ssize_t id; 
  LAVRow curBlock, prevBlock;

  ssize_t refGapLen, qryGapLen;
  id = 0;
  ssize_t curRow = 0;
  // Process every row
  // Get first coordinates
  MYSQL_ROW row;
  row = mysql_fetch_row(result);
  GetBlock(row, curBlock);
  LAVRow oend, obegin, onone;
  onone.tStart = -1;
  onone.tEnd   = -1;
  onone.qStart = -1;
  onone.qEnd   = -1;
  while (row = mysql_fetch_row(result)) {
    //    std::cout << "got row: " << row[0] << std::endl;
    prevBlock = curBlock;
    GetBlock(row, curBlock);

    // Look for a gap in the two blocks.
    // The gap has to be greater than min and less than max.

    // simple check for now.  Expand to something that allows for
    // a gap in both sequences.
    refGapLen = curBlock.tStart - prevBlock.tEnd + 1;
    qryGapLen = curBlock.qStart - prevBlock.qEnd + 1;
    std::cout << "refgaplen: " << refGapLen << " " << qryGapLen << std::endl;
    if (refGapLen > minGapLength &&
	refGapLen < maxGapLength) {
      std::cout << "found a gap at " << curBlock.tStart 
		<< " - " << prevBlock.tEnd + 1 
		<< "   " << curBlock.qStart << " - " << prevBlock.qEnd + 1 << std::endl;

      // Found a gap.  Now look for it in the other alignments.
      ssize_t s;
      // Positions in the species, and species query.
      ssize_t sStart, sEnd, sqStart, sqEnd, stStart, stEnd;
      for (s = 0; s < species.size(); s++) {
	std::cout <<"checking: " << species[s] << std::endl;
	tableName = species[s] + "_" + qrySpecies + "_" + seqName;
	
	// Find the block that contains the end of the previous block;
	obegin = onone;
	LookupReferencePosition(query, tableName, prevBlock.qEnd, stStart);
	std::cout << tableName << " lrp: " << prevBlock.qEnd << " " << stStart << std::endl;
	LookupReferencePosition(query, tableName, curBlock.qStart, stEnd);
	std::cout << tableName << " lrp: " << curBlock.qStart << " " << stEnd << std::endl;
	query << "SELECT * from " << tableName << " WHERE tStart <= " << stStart 
	      << " AND tEnd >= " << stStart;
	query.execute();
	if ((result = mysql_store_result(query.conn)) != NULL) {
	  numRows = mysql_num_rows(result);
	  std::cout << "got : " << numRows << std::endl;
	  if ( numRows > 0) {
	    row = mysql_fetch_row(result);
	    std::cout << "row: " << row[0] << " " << row[1] << " " << row[2] << " "
		      << row[3] << " " << row[4] << std::endl;
	    GetBlock(row, obegin);
	  }
	}	  
	
	oend = onone;
	query << "SELECT * from " << tableName << " WHERE tEnd <= " << stEnd
	      << " AND tStart >= " << stEnd << std::endl;
	if ((result = mysql_store_result(query.conn)) != NULL) {
	  numRows= mysql_num_rows(result);
	  std::cout << "got : " << numRows << std::endl;
	  if ( numRows > 0) {
	    row = mysql_fetch_row(result);
	    std::cout << "row: " << row[0] << " " << row[1] << " " << row[2] << " "
		      << row[3] << " " << row[4] << std::endl;
	    GetBlock(row, oend);
	  }
	}
	else {
	  oend = onone;
	}
	
	if (obegin.tStart != -1) {
	  std::cout << "got beginning block: " << obegin.tStart << " " << obegin.tEnd 
		    << " " << obegin.tStart << " " << obegin.tEnd << std::endl;
	}
	if (oend.tStart != -1) {
	  std::cout << "got beginning block: " << oend.tStart << " " << oend.tEnd 
		    << " " << oend.qStart << " " << oend.qEnd << std::endl;
	}
      }
    }
  }
}


void InitEnv(int argc, char* argv[], 
	     std::string &dbName,
	     std::string &refSpecies,
	     std::string &qrySpecies,
	     std::string &sequence,
	     std::vector<std::string> &species,
	     ssize_t &minGapLength,
	     ssize_t &maxGapLength,
	     ssize_t &minQryGapLength,
	     ssize_t &maxQryGapLength,
	     std::string &seqTableName) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "d:s:m:M:q:Q:")) != EOF){
    switch(copt) {
    case 's':
      seqTableName = optarg;
      continue;
    case 'd':
      dbName = optarg;
      continue;
    case 'q':
      minQryGapLength = atoi(optarg);
      continue;
    case 'Q':
      maxQryGapLength = atoi(optarg);
      continue;
    case 'm':
      minGapLength = atoi(optarg);
      continue;
    case 'M':
      maxGapLength = atoi(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  ssize_t i = optind;
  if (argc - i < 3) {
    std::cout <<" not enough arguments " << optind << std::endl;
    PrintUsage();
    exit(0);
  }
  refSpecies = argv[i];
  i++;
  qrySpecies = argv[i];
  i++;
  sequence = argv[i];
  i++;
  while (i < argc) {
    species.push_back(argv[i]);
    i++;
  }
}

void PrintUsage() {
  std::cout << "fog: program to find orthologous gaps given a bunch of pairwise alingments " << std::endl;
  std::cout << "usage: fog [-d database(lav_alignments)] [-s seqTableName(sequences)] [-m min] [-M max] refSpecies querySpeceis sequence species2 species3... "
	    << std::endl;
  std::cout << "-m minGapLength " << std::endl;
  std::cout << "-M maxGapLength " << std::endl;
  std::cout << "-d database " << std::endl;
  std::cout << "-s seqTableName " << std::endl;
}
