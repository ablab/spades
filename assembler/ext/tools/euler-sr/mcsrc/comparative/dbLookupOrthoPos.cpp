/***************************************************************************
 * Title:          dbLookupOrthoPos.cpp 
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
#include <set>

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

ssize_t ReadDBDescription(std::string       dbDescriptionFile,
		      DBDescriptionMap &dbDescriptions);

void PrintUsage() {
  std::cout << "dbloop: database based lookup-orthologous pos.  Extrapolate" 
	    << std::endl;
  std::cout << " between gaps in the alignment. " << std::endl;	
  std::cout << "usage: dbloop  [-f config_file] [-n] [-d database -s seqNameTable] [-p posfile|-l fromSpec toSpec region fromPos]." 
	    << std::endl;
  std::cout << "the format of the posfile is: " << std::endl;
  std::cout << "target ref qry pos chrom dir " << std::endl;
  std::cout << "target is the species to find the postion in. " << std::endl;
  std::cout << "ref is the original species " << std::endl;
  std::cout << "qry is ignored.  chrom is the sequence to look in, dir should be 0 " << std::endl;
}

class LookupQuery {
public:
  std::string refSeq;
  std::string qrySeq;
  std::string region;
  ssize_t refPos;
  ssize_t strand;
};
  
ssize_t InitEnv(int argc, char* argv[], 
	    std::string &posFile,
	    std::string &dbName,
	    std::string &defaultsFile,
	    std::string &outFileName,
	    std::string &seqTableName,
	    ssize_t &skipRevLookup,
	    ssize_t &transRevStrand,
	    ssize_t &verbosity,
	    LookupQuery &query) {
  ssize_t copt;
  ssize_t i;
  ssize_t posSpecified = 0;
  verbosity = 0;
  while ( (copt=getopt(argc, argv, "f:p:o:d:s:l:t:nv")) != EOF){
    switch(copt) {
    case 'v':
      ++verbosity;
      continue;
    case 'n':
      skipRevLookup = 1;
      continue;
    case 'f':
      defaultsFile = optarg;
      continue;
    case 't':
      transRevStrand = 1;
      continue;
    case 'o':
      outFileName = optarg;
      continue;
    case 'p':
      posFile = optarg;
      posSpecified = 1;
      continue;
    case 'd':
      dbName = optarg;
      continue;
    case 's':
      seqTableName = optarg;
      continue;
    case 'l':
      query.refSeq = optarg;
      query.qrySeq = argv[optind++];
      query.region = argv[optind++];
      query.refPos = atoi(argv[optind++]);
      query.strand = atoi(argv[optind++]);
      posSpecified = 1;
      continue;
    default:
      PrintUsage();
      exit(1);
      continue;
    }
  }
  if (!posSpecified ) {
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
		std::vector<ssize_t> &positions,
		std::vector<std::string> &targetSpecies,
		std::vector<std::string> &refSpecies,
		std::vector<std::string> &qrySpecies,
		std::vector<std::string> &chromosomes,
		std::vector<ssize_t> &dir); 


int main(int argc, char* argv[]) {
  std::string posFileName;
  std::string outputFileName;
  std::string seqTableName;
  std::string dbName;
  std::string defaultsFile = "";
  ssize_t transRevStrand = 0;
  ssize_t skipRevLookup = 0;
  dbName       = "lav_alignments";
  seqTableName = "sequences";
  LookupQuery lookupQuery;
  ssize_t verbosity;
  InitEnv(argc, argv, posFileName, 
	  dbName,
	  defaultsFile,
	  outputFileName,
	  seqTableName, 
	  skipRevLookup, transRevStrand,
	  verbosity,
	  lookupQuery);

  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName, defaultsFile);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);
  std::vector<ssize_t> pos;
  std::vector<std::string> chrom;
  std::vector<std::string> targetSpecies;
  std::vector<std::string> refSpecies;
  std::vector<std::string> qrySpecies;
  std::vector<ssize_t> dir;
  if (posFileName != "") {
    ReadPosFile(posFileName, pos, targetSpecies, refSpecies, qrySpecies, chrom, dir);
  }
  else {
    targetSpecies.push_back(lookupQuery.qrySeq);
    refSpecies.push_back(lookupQuery.refSeq);
    chrom.push_back(lookupQuery.region);
    dir.push_back(lookupQuery.strand);
    pos.push_back(lookupQuery.refPos);
  }
  std::string tableName, tmpTableName, chainTableName, netTableName;

  ssize_t p;
  ssize_t newPos;
  ssize_t refStart, refEnd, qryStart, qryEnd, qryStrand, chainId;
  ssize_t qryLen;
  ssize_t refLen;
  ssize_t numRows;
  std::set<std::string> tmpTables;
  
  ssize_t rowFound;
  MYSQL_RES *res;
  MYSQL_ROW row;
  ssize_t refPos, refStrand;
  refStrand = 0;

  for (p = 0; p < pos.size(); p++) {
    // Get the length of the sequence
    //    std::cout << pos[p] << " " << refSpecies[p] << " " << chrom[p] << " " << dir[p] << std::endl;
    newPos = pos[p];
    char pidStr[100];
    tableName = targetSpecies[p] + "_" + refSpecies[p] + "_" + chrom[p];
    sprintf(pidStr, PRI_PID, getpid());
    tmpTableName = tableName + "_tmp_" + pidStr; 
    chainTableName = tableName + "_chain";
    netTableName   = tableName + "_net";
    newPos = pos[p];
    ssize_t retVal;
    ssize_t deleted;
    if (!GetSequenceLength(query, seqTableName, refSpecies[p], chrom[p], refLen)) {
      std::cout  << "could not get sequence length for " << refSpecies[p] << std::endl;
      exit(1);
    }
    refStrand = 0;
    ssize_t orthStrand = 0;
    ssize_t orthPos = -1;
    retVal  = LookupOrthologousPosition(query, 
					tableName, 
					refStrand, pos[p], refLen,
					orthPos, orthStrand, tmpTables, 1024, deleted, skipRevLookup, verbosity);

    if (retVal != 0) {
      std::cout << refSpecies[p] << " " << targetSpecies[p] << " " << orthPos << " " << orthStrand <<  std::endl;
    }
    else {
      std::cerr << refSpecies[p] << " pos: " << pos[p] << " not found " << std::endl;
    }
  }
}
  

ssize_t ReadPosFile(std::string posFileName, 
		std::vector<ssize_t> &positions,
		std::vector<std::string> &targetSpecies,
		std::vector<std::string> &refSpecies,
		std::vector<std::string> &qrySpecies,
		std::vector<std::string> &chromosomes,
		std::vector<ssize_t> &dirVect) {
  std::ifstream posFile;
  openck(posFileName, posFile);
  ssize_t pos, dir;
  std::string chrom, ref, qry, target;
  while (posFile) {
    if (! (posFile >> target >> ref >> qry >> pos >> chrom >> dir)) 
      break;
    positions.push_back(pos);
    targetSpecies.push_back(target);
    refSpecies.push_back(ref);
    qrySpecies.push_back(qry);
    chromosomes.push_back(chrom);
    dirVect.push_back(dir);
  }
  posFile.close();
}


ssize_t ReadDBDescription(std::string dbDescriptionFile,
		      DBDescriptionMap &dbDescriptions) {

  std::ifstream dbDesFile;
  
  openck(dbDescriptionFile, dbDesFile);
  std::string specName, dbName;
  std::string qryName, qryTableName;
  ssize_t nTab;
  ssize_t n;
  ssize_t badFile = 0;
  while (dbDesFile && ! badFile) {
    if (! (dbDesFile >> specName >> dbName >> nTab)) 
      break;
    DBDescription *dbDescription = new DBDescription;
    dbDescriptions[specName] = dbDescription;
    dbDescription->dbName = dbName;
    for (n = 0; n < nTab; n++) {
      if (! (dbDesFile >> qryName >> qryTableName)) {
	badFile = 1;
	break;
      }
      dbDescription->specTable[qryName] = qryTableName;
    }
  }
}
