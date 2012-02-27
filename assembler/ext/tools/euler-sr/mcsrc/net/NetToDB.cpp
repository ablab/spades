/***************************************************************************
 * Title:          NetToDB.cpp 
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
#include <sstream>
#include "net/NetReader.h"
#include "blocks/BlockDB.h"
#include "blocks/BlockLib.h"
#include "mysql/mysql.h"

int main(int argc, char *argv[]) {

  std::string netFileName;
  std::string dbName;
  std::string netTableName;
  std::string refSpecies, qrySpecies, sequence;
  std::string defaultsFile;
  if (argc < 6) {
    std::cout << "usage: netToDB netFileName dbName refSpecies qrySpecies sequence [ defaultsFile ]" 
	      << std::endl;
    exit(0);
  }
  netFileName = argv[1];
  dbName      = argv[2];
  refSpecies  = argv[3];
  qrySpecies  = argv[4];
  sequence    = argv[5];
  if (argc == 7) 
    defaultsFile = argv[6];

  std::stringstream netTableNameStrm;
  netTableNameStrm << refSpecies << "_" << qrySpecies << "_" << sequence << "_net";
  netTableName = netTableNameStrm.str();
  std::cout << "defaults file: " << defaultsFile << std::endl;
  MYSQL *conn;
  ConnectToDB(conn, dbName.c_str(), defaultsFile);

  DBQuery query;
  query.conn = conn;

  NetFile netFile;
  NetReader::ReadNetFile(netFileName, netFile);

  std::vector<std::string> columns, indices;
  columns.push_back("id"); columns.push_back(" ssize_t AUTO_INCREMENT PRIMARY KEY ");
  columns.push_back("level"); columns.push_back(" INT ");
  columns.push_back("score"); columns.push_back(" BIGINT ");
  columns.push_back("orientation"); columns.push_back("TINYINT");
  columns.push_back("class"); columns.push_back("TINYINT");
  columns.push_back("chainid"); columns.push_back(" INT ");
  columns.push_back("tStart"); columns.push_back(" INT ");
  columns.push_back("tSize"); columns.push_back(" INT ");
  columns.push_back("qStart"); columns.push_back(" INT ");
  columns.push_back("qSize"); columns.push_back(" INT ");
  indices.push_back("chainid");

  CreateTable(query, netTableName, columns, indices);

  ssize_t i;
  for (i = 0; i < netFile.size(); i++) {
    if (netFile.nets[i]->chainClass == Net::fill) {
      if (netFile.nets[i] < 0) {
	std::cout << "fill specified without id, quitting " << std::endl;
	exit(0);
      }
      query << "INSERT INTO " << netTableName 
	    << " (level, score, orientation, class, tStart, tSize, qStart, qSize, chainid) VALUES  (" 
	    << netFile.nets[i]->level << "," 
	    << netFile.nets[i]->score << "," 
	    << netFile.nets[i]->orientation << "," 
	    << netFile.nets[i]->chainClass << "," 
	    << netFile.nets[i]->tStart << "," 
	    << netFile.nets[i]->tSize << "," 
	    << netFile.nets[i]->qStart << "," 
	    << netFile.nets[i]->qSize << "," 
	    << netFile.nets[i]->id << ")";
      query.execute();
    }
  }
}
