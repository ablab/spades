/***************************************************************************
 * Title:          MakeSelfTable.cpp 
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
    std::cout << "usage: maketemp [-f options_file] dbname seqnem spec1 [spec2 ...] " << std::endl;
    exit(1);
  }
  std::string dbName, seqName, optionsFile;
  optionsFile = "";
  int argi = 0;
  if (strcmp(argv[1], "-f") == 0) {
    optionsFile = argv[2];
    argi = 3;
  }
  else {
    argi = 1;
  }
  dbName = argv[argi++];
  seqName = argv[argi++];
  
    // Connect to the database;
  MYSQL *dbConn;
  std::cout << "using options: " << optionsFile << std::endl;
  ConnectToDB(dbConn, dbName, optionsFile);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);


  std::string tableName, netTableName, chainTableName, tempTableName;
  std::string selfTable;
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
  for (s1 = 0; s1 < species.size(); s1++) {
    tableName = species[s1] + "_" + species[s1] + "_" + seqName;
    chainTableName = tableName + "_chain";
    netTableName   = tableName + "_net";
    selfTable = species[s1] + "_" + seqName + "_self";
    std::cout << "making table: " << selfTable << std::endl;
    query << "CREATE TABLE IF NOT EXISTS " << selfTable << " AS SELECT " 
	  << tableName << ".id, " 
	  << tableName << ".tStart, " 
	  << tableName << ".tEnd, "
	  << tableName << ".qStart, "
	  << tableName << ".qEnd, "
	  << tableName << ".strand," 
	  << chainTableName << ".chainId " 
	  << " FROM " << tableName << ", " << chainTableName 
	  << " WHERE " << tableName << ".id = " << chainTableName << ".lavid AND "
	  << tableName << ".strand = 1 ";
    query.execute();
  }
  
}
