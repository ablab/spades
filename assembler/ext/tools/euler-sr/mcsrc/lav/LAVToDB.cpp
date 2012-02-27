/***************************************************************************
 * Title:          LAVToDB.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std library
#include <iostream>
#include <fstream>
#include <string>
#include <exception>
#include <stdio.h>
#include <sstream>

// other library
#include "mysql/mysql.h"

// mine
#include "utils.h"
#include "lav/LAVTable.h"
#include "lav/LAVFile.h"
#include "lav/LAVAlignedContig.h"
#include "lav/LAVBlock.h"
#include "lav/LAVReader.h"
#include "lav/LAVUtils.h"

#include "blocks/BlockDB.h"
#include "blocks/BlockLib.h"


void CreateTableName(std::string &tName, 
		     std::string &qName, 
		     std::string &seqName, 
		     std::string &tableName) {
  ssize_t tIndex, qIndex;
  tableName = tName + "_" + qName + "_" + seqName;
}


int main(int argc, char* argv[]) {
  std::string lavFileName, refSpecies, qrySpecies, seqName, database, defaultsFile;
  std::string tabFileName;
  std::string seqTable;
  ssize_t qrySeqId, refSeqId;
  ssize_t useDB = 1;
  defaultsFile = "";
  refSpecies = "";
  qrySpecies = "";
  seqName    = "";
  database   = "alignments";
  
  if (argc < 4) {
    std::cout << "usage: lavToDB [-f defaultsFile] lavFile refSepecies qrySpecies sequence -t |  database seqtable" << std::endl;
    return 1;
  }
  tabFileName = "";
  int argi;
  if (strcmp(argv[1], "-f") == 0) {
    defaultsFile = argv[2];
    argi = 3;
  }
  else {
   argi = 1;  
  }
  lavFileName = argv[argi++];
  refSpecies  = argv[argi++];
  qrySpecies  = argv[argi++];
  seqName     = argv[argi++];

  if (strcmp(argv[argi], "-t") == 0) {
    useDB = 0;
  }
  else {
    database    = argv[argi++];
    seqTable    = argv[argi++];
  }

  std::ofstream tabOutFile;
  // Connect to the database
  MYSQL *dbConn;
  if (useDB) {
    ConnectToDB(dbConn, database, defaultsFile);
  }

  //  ParseFileName(lavFileName, refSpecies, qrySpecies, seqName);

  
  DBQuery query(dbConn);

  // Get the information to add to the table.
  // the ids of the sequences the alignment corresponds to
 
  refSeqId = 0;
  qrySeqId = 1;
  
  std::string tableName;
  CreateTableName(refSpecies, qrySpecies, seqName, tableName);
  
  tabFileName = tableName + ".txt";
  if (!useDB) {
    std::cout << "opening " << tabFileName << std::endl;
    openck(tabFileName, tabOutFile, std::ios::out);
  }

  LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

  ssize_t a, b, cb;
  LAVAlignedContig *alignedContig;
  double score;
  std::string tName, qName;
  ssize_t tSize, tStart, tEnd, qSize, qStart, qEnd;
  char tStrand, qStrand;
  ssize_t size, dt, dq;
  std::string word;
  ssize_t good = 1;
  ssize_t chainStarted = 0;
  ssize_t tPos, qPos;
  ssize_t dbInitialized = 0;

  std::vector<std::string > columns, indices;
  columns.push_back("id"); columns.push_back(" ssize_t AUTO_INCREMENT PRIMARY KEY ");
  columns.push_back("tStart"); columns.push_back("ssize_t");     
  columns.push_back("tEnd");   columns.push_back("ssize_t");     
  columns.push_back("qStart"); columns.push_back("ssize_t");     
  columns.push_back("qEnd");   columns.push_back("ssize_t");     
  columns.push_back("strand"), columns.push_back("ssize_t");
  columns.push_back("chainId"); columns.push_back("ssize_t");     
  columns.push_back("tid"); columns.push_back("ssize_t");
  columns.push_back("qid"); columns.push_back("ssize_t");
  indices.push_back("tStart");
  indices.push_back("tEnd");
  indices.push_back("qStart");
  indices.push_back("qEnd");
  indices.push_back("chainId");


  ssize_t strand;
  if (useDB) {
    CreateTable(query, tableName, columns, indices);
    query << "LOCK TABLES " << tableName << " WRITE ";
    query.execute();
  }

  ssize_t chainId = 1;
  for (a = 0; a < lavFile.alignments.size(); a++) {
    //    std::cout << a << " /// " << lavFile.size() << std::endl;
    alignedContig = lavFile.alignments[a];
    LAVBlock *block;
    strand = alignedContig->qryContig.strand;
    for (b = 0; b < alignedContig->alignments.size(); b++) {
      block = alignedContig->alignments[b];
      //	std::cout << b << " / " << alignedContig->alignments.size() << std::endl;
      if (useDB) {
	for (cb = 0; cb < block->size(); cb++) { 
	  ++chainId;
	  if (useDB) {
	    EnqueueStoreRow(query, tableName, 
			    block->refALBegin[cb],
			    block->refALEnd[cb],
			    block->qryALBegin[cb],
			    block->qryALEnd[cb], 
			    strand,
			    chainId,
			    refSeqId, qrySeqId, 1000);
	  }
	}
      }
      else {
	block->TabbedPrintBlock(tabOutFile, chainId, strand, refSeqId, qrySeqId);
      }
    }
    if (useDB) {
      query.execute();
    }
  }
  if (useDB) {
    query << "UNLOCK TABLES ";
    query.execute();
  }
  else {
    tabOutFile.close();
  }
}
 
