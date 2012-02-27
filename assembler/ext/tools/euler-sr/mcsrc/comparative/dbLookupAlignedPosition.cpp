/***************************************************************************
 * Title:          dbLookupAlignedPosition.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/25/2009
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
#include <map>

// 3rd party includes
#include "mysql/mysql.h"

// my stuff
#include "utils.h"
#include "lav/LAVTable.h"
#include "blocks/BlockDB.h"
#include "blocks/dbOrthoPos.h"

#include "compatibility.h"
	     
class DBDescription {
public:
  std::string dbName;
  std::map<std::string, std::string> specTable;
};

typedef std::map<std::string, DBDescription*> DBDescriptionMap;


void PrintUsage() {
  std::cout << "dbloop: database based lookup pos.  Do not extrapolate " 
	    << std::endl;
  std::cout << "between gaps in the alignment. " << std::endl;
  std::cout << "usage: dblap  [-d database -s seqNameTable] posfile." << std::endl;
  std::cout << "posfile contains a line with each quer species to consider, then a " << std::endl
	    << "line with each target species coordinate, and direction. " << std::endl
	    << "example:"  << std::endl
	    << "mouse rat platypus " << std::endl
	    << "human ENr222  100000 0" << std::endl
	    << "dog  ENm011 123456 0" << std::endl;
}

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &posFile,
	    std::string &dbName,
	    std::string &outFileName,
	    std::string &seqTableName,
	    ssize_t &transRevStrand,
	    ssize_t &verbosity) {
  ssize_t copt;
  ssize_t i;
  verbosity =0;
  while ( (copt=getopt(argc, argv, "o:d:s:tv")) != EOF){
    switch(copt) {
    case 'v':
      verbosity++;
      continue;
    case 't':
      transRevStrand = 1;
      continue;
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
    posFile = argv[i++];
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

ssize_t ReadPosFile(std::string posFileName, 
		std::vector<std::string> &querySpecies,
		std::vector<std::string> &refSpecies,
		std::vector<std::string> &refSeq,
		std::vector<ssize_t> &refPositions);

int main(int argc, char* argv[]) {
  std::string posFileName;
  std::string outputFileName;
  std::string seqTableName;
  std::string dbName;
  ssize_t transRevStrand = 0;
  dbName       = "lav_alignments";
  seqTableName = "sequences";
  ssize_t verbosity = 0;
  InitEnv(argc, argv, posFileName, 
	  dbName,
	  outputFileName,
	  seqTableName, transRevStrand, verbosity);

  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);

  std::vector<ssize_t> pos;
  std::vector<std::string> chrom;
  std::vector<std::string> targetSpecies;
  std::vector<std::string> refSpecies;
  std::vector<std::string> qrySpecies;
  std::vector<ssize_t> dir;
  ReadPosFile(posFileName, qrySpecies, refSpecies, chrom, pos);
  std::string tableName, tmpTableName, chainTableName, netTableName;

  ssize_t p, q;
  ssize_t newPos;
  ssize_t refStart, refEnd, qryStart, qryEnd, qryStrand, chainId;
  ssize_t qryLen;
  ssize_t refLen;
  ssize_t numRows;
  std::set<std::string> tmpTables;
  
  for (p = 0; p < pos.size(); p++) {
    // Get the length of the sequence
    query << "SELECT length from " << seqTableName << " where name =" << dbstr(refSpecies[p]) 
	  << " AND sequence="<< dbstr(chrom[p]);
    GetOneInt(query, refLen);


    ssize_t rowFound;
    MYSQL_RES *res;
    MYSQL_ROW row;
    ssize_t refPos;
    char pidStr[100];
    ssize_t refStrand;
    ssize_t deleted;
    for (q = 0; q < qrySpecies.size(); q++ ) {
      tableName = refSpecies[p] + "_" + qrySpecies[q] + "_" + chrom[p];
      sprintf(pidStr, PRI_PID, getpid());
      tmpTableName = tableName + "_tmp_" + pidStr; 
      chainTableName = tableName + "_chain";
      netTableName   = tableName + "_net";
      newPos = pos[p];


      rowFound = LookupOrthologousPosition(query, 
					   tableName, 
					   0, pos[p], refLen, refPos, refStrand,
					   tmpTables, 100, deleted, 0, verbosity);

      if (deleted) {	
	// Do we care for now that the position was deleted????
	// No.
      }
    } // end looking in + strand
  } // end looking through each species
}
  

ssize_t ReadPosFile(std::string posFileName, 
		std::vector<std::string> &querySpecies,
		std::vector<std::string> &refSpecies,
		std::vector<std::string> &refSeq,
		std::vector<ssize_t> &refPositions) {
  std::ifstream posFile;
  openck(posFileName, posFile);
  std::string spec;
  ssize_t pos;
  std::string seq, refSpec;
  // read the species to lookup
  char line[1000];
  posFile.getline(line, 1000);
  char *tok;
  tok = strtok(line, " \t");
  while ( tok != NULL ) {
    querySpecies.push_back(tok);
    std::cout << "got qry species: " << tok << std::endl;
    tok = strtok(NULL, " \t\n");
  }
  
  while (posFile) {
    if (! (posFile >> refSpec >> seq >>  pos)) 
      break;
    refSpecies.push_back(refSpec);
    refSeq.push_back(seq);
    refPositions.push_back(pos);
  }
  posFile.close();
  std::cout << "done reading positions " << std::endl;
}


