/***************************************************************************
 * Title:          dbInvCheck.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
// std includes
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <unistd.h>

// 3rd party includes
#include "mysql/mysql.h"

// my stuff
#include "utils.h"
#include "lav/LAVTable.h"
#include "blocks/BlockDB.h"

#include "compcommon.h"
	     
class DBDescription {
public:
  std::string dbName;
  std::map<std::string, std::string> specTable;
};

typedef std::map<std::string, DBDescription*> DBDescriptionMap;


void PrintUsage() {
  std::cout << "dbloop: database based lookup-orthologous pos." << std::endl;
  std::cout << "usage: dbloop  [-d database -s seqNameTable] posfile." << std::endl;
}

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &dbName,
	    std::string &lavFile,
	    std::string &outFileName,
	    std::string &seqTableName) {
  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "o:d:s:")) != EOF){
    switch(copt) {
    case 'o':
      outFileName = optarg;
      continue;
    case 'd':
      dbName = optarg;
      continue;
    case 's':
      seqTableName = optarg;
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


int main(int argc, char* argv[]) {
  std::string posFileName;
  std::string outputFileName;
  std::string seqTableName;
  std::string dbName;
  std::string lavFile;
 
  std::string refSpec, qrySpec, chrom;
  dbName       = "lav_alignments";
  seqTableName = "sequences";
  InitEnv(argc, argv, posFileName, 
	  dbName,
	  lavFile,
	  outputFileName,
	  seqTableName);

  ParseFileName(lavFileName, refSpec, qrySpec, chrom);
  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);
  std::vector<ssize_t> pos;
  
  std::string tableName;
  
  std::vector<ssize_t> invStart, invEnd;

  
  ssize_t p;
  ssize_t newPos;
  ssize_t refStart, refEnd, qryStart, qryEnd, qryStrand, chainId;
  ssize_t qryLen;

  
  
  for (p = 0; p < pos.size(); p++) {
    // Get the length of the sequence

    query << "SELECT length from " << seqTableName;
    GetOneInt(query, qryLen);
    newPos = qryLen - pos[p] + 1; // don't index directly inside a repeat.
    tableName = refSpecies[p] + "_" + qrySpecies[p] + "_" + chrom[p];

    std::cout << "selecting from table: " << tableName << std::endl;
    query << "select tStart, tEnd, qStart,qEnd,strand,chainID from "
	  << tableName << " where qStart < " << newPos << " and qEnd > " 
	  << newPos;
    std::cout << "selecting: " << query.preview() << std::endl;
    if (query.execute()) {
      HandleError(query.conn);
    }
    MYSQL_RES *res;
    MYSQL_ROW row;
    if ((res = mysql_store_result(query.conn)) == NULL)
      HandleError(query.conn);
    
    ssize_t numRows = mysql_num_rows(res);
    ssize_t refPos;
    ssize_t rowFound = 0;
    while (numRows > 0 &&
	   FetchRow(res,refStart, refEnd, qryStart, qryEnd, 
		    qryStrand, chainId)) {
      rowFound = 1;
      GetOffsetPosition(qryStart, qryEnd, pos[p] +1,
			refStart, refEnd, refPos);
      
      std::cout << refSpecies[p] << " " << qrySpecies[p] 
		<< pos[p] + 1 << " " << refPos << std::endl;
      numRows--;
    }
    if (rowFound == 0) {
      std::cout << "Didn't find a match for " << newPos << std::endl;
    }
  }
}

