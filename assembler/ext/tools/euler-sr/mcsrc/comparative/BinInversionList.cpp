/***************************************************************************
 * Title:          BinInversionList.cpp 
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
#include <vector>
#include <set>
#include <map>
#include <sstream>

#include "mysql/mysql.h"

#include "utils.h"
#include "blocks/BlockDB.h"
#include "net/NetDB.h"
#include "blocks/dbOrthoPos.h"
#include "CharTree.h"
#include "InversionDB.h"
#include "InversionBins.h"
#include "InversionUtils.h"

int main(int argc, char* argv[]) {

  std::string dbName, validInvFileName, binFileName, seqName;
  std::string graphName;
  validInvFileName = "";
  binFileName = "";
  std::cout << "THIS PROGRAM HAS BEEN MOTHBALLED" << std::endl;
  exit(0);
  if (argc != 5) {
    std::cout << "usage: bininversions dbName sequenceName validatedInvFile binFileName " << std::endl;
    exit(1);
  }
  dbName           = argv[1];
  seqName          = argv[2];
  validInvFileName = argv[3];
  binFileName      = argv[4];

  
  // Connect to the database 
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);

  StringVector species;
  InversionList invList;
  BinMap binnedInversions;
  std::vector<ssize_t> startPos, endPos;
  ssize_t deleted;
  // Read in the inversions and where they belong.
  ReadValidatedFile(validInvFileName, species, invList);
  BinInversionList(query, 
		   invList, 
		   seqName,
		   startPos, endPos,
		   binnedInversions);

  std::ofstream outfile;
  openck(binFileName, outfile);
  PrintBins(species, binnedInversions, outfile);

  return 0;
}
