/***************************************************************************
 * Title:          FixSeqLengths.cpp 
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

#include "mysql/mysql.h"

#include "DNASequence.h"
#include "SeqReader.h"
#include "blocks/BlockDB.h"

void ParseFileName(std::string fileName, 
		   std::string &refName, 
		   std::string &seqName) {
  // Find the ending component of the path.
  ssize_t startRefInd, endRefInd, qryStartInd, qryEndInd;

  startRefInd = fileName.rfind("/");
  if (startRefInd < 0) 
    startRefInd = 0;
  
  endRefInd = fileName.find(".EN", startRefInd);
  if (endRefInd == fileName.npos) {
    std::cout << "error parsing " << fileName << std::endl;
    exit(0);
  }
  refName = fileName.substr(startRefInd, endRefInd - startRefInd );
  ssize_t faStart;
  faStart = fileName.find(".fa", endRefInd);
  if (faStart == fileName.npos) {
    std::cout << "error parsing " << fileName << std::endl;
    exit(0);
  }
  seqName = fileName.substr(endRefInd + 1, faStart  - endRefInd -1);
  std::cout << "got rn: " << refName << " sn: " << seqName << std::endl;
}
  


int main(int argc, char *argv[]) {

  MYSQL *conn;
  
  std::string database, fileName;
  std::string tableName, alignTableName;
  if (argc != 4) {
    std::cout << "usage: fixseqlengths database tablename filename " << std::endl; 
    exit(1);
  }
  database = argv[1];
  tableName= argv[2];
  fileName = argv[3];
  std::string refName, qryName, seqName;
  ParseFileName(fileName, refName, seqName);
  std::vector<std::string> columns, indices;
  
  columns.push_back("id"); columns.push_back("INT AUTO_INCREMENT PRIMARY KEY");
  columns.push_back("name"); columns.push_back("varchar(255)");
  columns.push_back("sequence"); columns.push_back("varchar(255)");
  columns.push_back("file"); columns.push_back("varchar(255)");
  columns.push_back("length"); columns.push_back("ssize_t");

  ConnectToDB(conn, database.c_str());

  DBQuery query(conn);
  CreateTable(query, tableName, columns, indices, 0);

  DNASequence seq;
  SeqReader::GetSeq(fileName, seq);
  query << "INSERT INTO " << tableName << " (name, sequence, file, length)" 
	<< "VALUES ( " 
	<< dbstr(refName) <<","<< dbstr(seqName) <<","<< dbstr(fileName) << "," << seq.length 
	<< ")";
  if (query.execute()) {
    std::cout << "error executing " << query.prevQuery << " " << mysql_error(query.conn) << std::endl;
    exit(0);
  }

  query << "SELECT MAX(id) from " << tableName;
  query.execute();

  MYSQL_RES *res;
  if ((res = mysql_store_result(query.conn)) == NULL) HandleError(conn);
  if (mysql_num_rows(res) == 0) {
    std::cout << "there must be at least a amx id " << std::endl;
    exit(0);
  }
  MYSQL_ROW row;
  if ((row = mysql_fetch_row(res)) == NULL) {
    std::cout <<"error: " << mysql_error(query.conn) << std::endl;
    exit(1);
  }
}
