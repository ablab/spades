/***************************************************************************
 * Title:          FixTableIndex.cpp 
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
#include "blocks/BlockDB.h";


int main(int argc, char *argv[]) {

  MYSQL *conn;
  
  std::string database, alignFileName;

  std::string alignmentTable, seqTable;
  if (argc != 4) {
    std::cout << "usage: fti database seqtable alignFileName " << std::endl; 
    exit(1);
  }
  database = argv[1];
  seqTable   = argv[2];
  alignFileName = argv[3];
  std::string refName, qryName, seqName;
  ParseFileName(alignFileName, refName, qryName, seqName);

  alignmentTable = refName + "_" + qryName + "_" + seqName;
  std::vector<std::string> columns, indices;

  ConnectToDB(conn, database.c_str());

  DBQuery query(conn);

  
  // Get the information to add to the table.
  // the ids of the sequences the alignment corresponds to
  query << "SELECT id from " << seqTable << " where name=" << dbstr(refName) << " and sequence="<< dbstr(seqName);
  ssize_t refSeqId, qrySeqId;
  GetOneInt(query, refSeqId);

  query << "SELECT id from " << seqTable << " WHERE name=" << dbstr(qryName) << " and sequence=" << dbstr(seqName);
  GetOneInt(query, qrySeqId);

  // update the alignment table
  query << "ALTER TABLE " << alignmentTable << " ADD COLUMN (tid ssize_t, qid ssize_t)";
  std::cout << "query executing: " << query.str() << std::endl;
  query.execute();
  /*
  if (query.execute()){ 
    std::cout << "Could not alter table " << alignmentTable << " " << mysql_error(query.conn);
    exit(1);
  }
  */
  // Update the values of the rows.
  ssize_t numRows;
  GetNumBlocks(query, alignmentTable, numRows);
  std::cout << "about to update: " << numRows << " rows " << std::endl;
  ssize_t i;
  query << "UPDATE " << alignmentTable << " SET tid=" << refSeqId << ", qid=" << qrySeqId << " WHERE id >= 0 and id <= " << numRows;
  query.execute();
}
