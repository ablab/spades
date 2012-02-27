/***************************************************************************
 * Title:          MakeTemp.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <glob.h>

#include "mysql/mysql.h"
#include "blocks/BlockDB.h"
#include "ValidatedInversion.h"
#include "InversionUtils.h"

int main(int argc, char* argv[]) {
  if (argc < 4) {
    std::cout << "usage: maketemp [-f options-file] dbname seqnem spec1 [spec2 ...] " << std::endl;
    exit(1);
  }
  std::string dbName, seqName, optionsFile;
  ssize_t verbosity= 0;
  ssize_t symmetric = 0;
  optionsFile = "";
  int argi = 1;
  while (strlen(argv[argi])> 0 and
	 argv[argi][0] == '-') {
    if (strcmp(argv[argi], "-f") == 0) {
      ++argi;
      optionsFile = argv[argi];
    }
    if (strcmp(argv[argi], "-v") == 0) {
      ++verbosity;
    }
    if (strcmp(argv[argi], "-S") == 0) {
      symmetric = 0;
    }
    ++argi;
  }
  dbName = argv[argi++];
  seqName = argv[argi++];
  
    // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName, optionsFile);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);


  std::string tableName, netTableName, chainTableName, tempTableName;
  std::string permTempTableName;

  StringVector species;
  StringSet  tempTableNames;
  std::string tempspec;
  std::string filename;
  glob_t g;
  while (argi < argc) {
    tempspec = argv[argi];
    species.push_back(tempspec);
    argi++;
  }
  ssize_t s1, s2;
  ssize_t lower, upper;
  if (symmetric)
    upper = species.size();
  else
    upper= species.size() - 1;

  for (s1 = 0; s1 < upper; s1++) {
    if  (symmetric) {
      lower = 0;
    }
    else {
      lower = s1+ 1;
    }
    for (s2 = lower; s2 < species.size(); s2++) {
      tableName = species[s1] + "_" + species[s2] + "_" + seqName;
      tempTableName  = tableName + "_temp1";
      chainTableName = tableName + "_chain";
      netTableName   = tableName + "_net";
      
      permTempTableName = tableName + "_temp";

      std::cout << " creating temporary table with: " << tableName 
		<< " " << chainTableName << " " << netTableName
		<< " " << tempTableName << std::endl;

      std::vector<ssize_t> firstLevel;
      std::string firstLevelP;
      GetFirstLevelNets(query, netTableName, firstLevel);
      BuildChainIdPredicate(chainTableName, firstLevel, firstLevelP);
      query << "DROP TABLE IF EXISTS " << permTempTableName;
      query.execute();

      query << "CREATE TABLE IF NOT EXISTS " << permTempTableName 
	    << " (id INTEGER, tStart INTEGER, tEnd INTEGER, " 
	    << " qStart INTEGER, qEnd INTEGER, chainId INTEGER , " 
	    << " strand INTEGER , KEY(tStart), KEY(tEnd), KEY(qStart), KEY(qEnd) )";
      query.execute();

      query << "INSERT INTO " << permTempTableName 
	    << " (id, tStart, tEnd, qStart, qEnd, chainId, strand) SELECT "
	    << tableName << ".id, " <<   tableName << ".tStart,"
	    << tableName << ".tEnd," <<  tableName << ".qStart, "
	    << tableName << ".qEnd," << tableName << ".chainId , " 
	    << tableName << ".strand FROM " << tableName << ", " 
	    << chainTableName << " WHERE " << tableName << ".chainId = " 
	    << chainTableName << ".lavId AND " << tableName << ".strand = 0 " 
	    << " AND " << firstLevelP;
      query.execute();


      // Now do the same for the reverse strand.

      std::string reverseName;
      reverseName = permTempTableName + "_reverse";
      GetFirstLevelNets(query, netTableName, firstLevel, 3, 1);
      BuildChainIdPredicate(chainTableName, firstLevel, firstLevelP);
      query << "DROP TABLE IF EXISTS " << reverseName;
      query.execute();

      query << "CREATE TABLE IF NOT EXISTS " << reverseName 
	    << " (id INTEGER, tStart INTEGER, tEnd INTEGER, " 
	    << " qStart INTEGER, qEnd INTEGER, chainId INTEGER , " 
	    << " strand INTEGER , KEY(tStart), KEY(tEnd), KEY(qStart), KEY(qEnd) )";
      query.execute();

      query << "INSERT INTO " << reverseName << " (id, tStart, tEnd, qStart, qEnd, chainId, strand) SELECT "
	    << tableName << ".id, " <<   tableName << ".tStart,"
	    << tableName << ".tEnd," <<  tableName << ".qStart, "
	    << tableName << ".qEnd," << tableName << ".chainId , " 
	    << tableName << ".strand FROM " << tableName << ", " 
	    << chainTableName << " WHERE " << tableName << ".chainId = " 
	    << chainTableName << ".lavId AND " << tableName << ".strand = 0 " 
	    << " AND " << firstLevelP;
      query.execute();
      
      
      std::cout << "used " << query.prevQuery << std::endl;
    }
  }
}
