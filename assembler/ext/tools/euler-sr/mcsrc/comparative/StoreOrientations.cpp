/***************************************************************************
 * Title:          StoreOrientations.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "mysql/mysql.h"
#include "blocks/BlockDB.h"



int main(int argc, char* argv[]) {

  if (argc <= 4) {
    std::cout << "usage: storient refspec qryspec seq " << std::endl;
    exit(0);
  }
  std::string dbName = "lav_alignments";
  std::string oTabName = "orientations";

  if (argc >= 5)
    dbName = argv[4];

  if (argc >= 6)
    oTabName = argv[5];

  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);

  
  std::vector<std::string> columns, indices;
  columns.push_back("id"); columns.push_back(" ssize_t AUTO_INCREMENT PRIMARY KEY ");
  columns.push_back("align_name"); columns.push_back("varchar(255)");
  columns.push_back("orientation"); columns.push_back("ssize_t");
  
  CreateTable(query, oTabName,columns, indices, 0);
  
  std::string refSpec, qrySpec, seq;
  refSpec = argv[1];
  qrySpec = argv[2];
  seq  = argv[3];
  std::string lavTableName = refSpec + "_" + qrySpec + "_" + seq;
  
  ssize_t forSize, revSize;
  query << "select sum (tEnd - tStart) from " << lavTableName << " where strand=0";
  
  if (!GetOneInt(query, forSize))
    forSize = -1;

  query << "select sum (tEnd - tStart) from " << lavTableName << " where strand=1";
  
  if (!GetOneInt(query, revSize))
    revSize = -1;
  
  ssize_t orientation;
  if (forSize > revSize) 
    orientation = 0;
  else
    orientation = 1;

  query << "INSERT INTO TABLE " << oTabName << " VALUES (0, " << dbstr(lavTableName) 
	<< ", " << orientation << ")";
  return query.execute();
} 
