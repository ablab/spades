/***************************************************************************
 * Title:          SelectNettedAlignments.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/25/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <unistd.h>
#include "mysql/mysql.h"

#include "compatibility.h"
#include "utils.h"
#include "blocks/BlockDB.h"
#include "lav/LAVReader.h"
#include "lav/LAVFile.h"
#include "lav/LAVUtils.h"
#include "lav/LAVPrinter.h"
#include "net/NetDB.h"

#include "compatibility.h"

int main(int argc, char *argv[]) {

  std::string lavFileName, refSpecies, qrySpecies, seqName, database;
  std::string outFileName;
  std::string seqTable;
  ssize_t qrySeqId, refSeqId;
  refSpecies = "";
  qrySpecies = "";
  seqName    = "";
  database   = "alignments";

  if (argc < 4) { 
    std::cout << "usage: selectNet database lavFile outfile [level]" << std::endl;
    exit(1);
  }
  database    = argv[1];
  lavFileName = argv[2];
  outFileName = argv[3];
  ssize_t level = 1;
  if (argc == 5) 
    level = atoi(argv[4]);

  std::string lavTableName, chainTableName, netTableName, tempTableName;

  // Connect to the database
  MYSQL *dbConn;

  ConnectToDB(dbConn, database);

  ParseFileName(lavFileName, refSpecies, qrySpecies, seqName);
  
  lavTableName = refSpecies + "_" + qrySpecies + "_" + seqName;
  chainTableName = lavTableName + "_chain";
  netTableName  = lavTableName + "_net";

  char pidStr[100];
  sprintf(pidStr, PRI_PID, getpid());
  tempTableName = lavTableName + "_tmp_" + pidStr;
  
  DBQuery query(dbConn);

  std::vector<ssize_t> firstLevel;
  std::string firstLevelP;
  GetFirstLevelNets(query, netTableName, firstLevel, level);
  BuildChainIdPredicate(chainTableName, firstLevel, firstLevelP);

  query << "CREATE TEMPORARY TABLE " << tempTableName << " AS SELECT "
	<< lavTableName << ".id FROM " << lavTableName << ", " << chainTableName 
	<< " WHERE " << lavTableName << ".id = " << chainTableName << ".lavId AND " 
	<< firstLevelP;

  query.execute();
  query << "ALTER TABLE " << tempTableName << " ADD INDEX (id) ";
  query.execute();
  std::cout << "done creating temporary table " << std::endl;
  ssize_t ac, a;
  LAVAlignedContig *alignedContig;
  ssize_t id = 0;
  ssize_t fetchedId;

  LAVFile lavFile;
  LAVReader::ReadLAVFile(lavFileName, lavFile);

  std::vector<LAVBlock*>::iterator acit, acend;
  std::vector<ssize_t>::iterator refBeginIt, refBeginEnd, refEndIt, 
    refEndEnd, qryBeginIt, qryBeginEnd, qryEndIt, qryEndEnd;
  ssize_t tStart, tEnd, qStart, qEnd;
  ssize_t idnum = 0;


  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    for (acit = alignedContig->alignments.begin();
	 acit != alignedContig->alignments.end(); ++acit) {
      refBeginIt = (*acit)->refALBegin.begin();
      refEndIt   = (*acit)->refALEnd.begin();
      qryBeginIt = (*acit)->qryALBegin.begin();
      qryEndIt   = (*acit)->qryALEnd.begin();
      
      for (; refBeginIt != (*acit)->refALBegin.end(); ++refBeginIt) {
	++idnum;
      }
    }
  }
  std::cout << "got " << idnum << " ids " << std::endl;
  for (ac = 0; ac < lavFile.alignments.size(); ac++) {
    alignedContig = lavFile.alignments[ac];
    for (acit = alignedContig->alignments.begin();
	 acit != alignedContig->alignments.end(); ) {
      refBeginIt = (*acit)->refALBegin.begin();
      refEndIt   = (*acit)->refALEnd.begin();
      qryBeginIt = (*acit)->qryALBegin.begin();
      qryEndIt   = (*acit)->qryALEnd.begin();
      
      for (; refBeginIt != (*acit)->refALBegin.end(); ) {
	++idnum;
	if (idnum % 1000 == 0)
	  std::cout << "processing " << idnum << std::endl;
	tStart = (*refBeginIt);
	tEnd   = (*refEndIt);
	qStart = (*qryBeginIt);
	qEnd   = (*qryEndIt);

	query << "SELECT id FROM " << lavTableName << " WHERE "
	      << "tStart = " << tStart << " AND " 
	      << "tEnd   = " << tEnd   << " AND "
	      << "qStart = " << qStart << " AND "
	      << "qEnd   = " << qEnd;
	GetOneInt(query, id);
	query << "SELECT id FROM " << tempTableName << " WHERE id = " << id;
	if (GetUpToInt(query, fetchedId) == 0) {
	  // Found an id that isn't on the main chain, get rid of it.
	  refBeginIt = (*acit)->refALBegin.erase(refBeginIt);
	  refEndIt   = (*acit)->refALEnd.erase(refEndIt);
	  qryBeginIt = (*acit)->qryALBegin.erase(qryBeginIt);
	  qryEndIt   = (*acit)->qryALEnd.erase(qryEndIt);
	}
	else {
	  ++refBeginIt;
	  ++refEndIt;
	  ++qryBeginIt;
	  ++qryEndIt;
	}
      }
      if ((*acit)->refALBegin.size() == 0) {
	acit = alignedContig->alignments.erase(acit);
      }
      else {
	++acit;
      }
    }
  }
  std::ofstream outfile;
  openck(outFileName, outfile, std::ios::out);
  LAVPrinter::PrintLAVFile(lavFile, outfile);
  return 0;
}
