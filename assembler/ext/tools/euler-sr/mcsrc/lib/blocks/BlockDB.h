/***************************************************************************
 * Title:          BlockDB.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BLOCK_DB_H_
#define BLOCK_DB_H_
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "mysql/mysql.h"

// my own stuff
#include "lav/LAVTable.h"

class DBQuery : public std::stringstream {
public:
  MYSQL *conn;
  std::string queryStr;
  std::string prevQuery;
  ssize_t queueSize;
  DBQuery() {
    queryStr = "";
    str(queryStr);
    conn = NULL;
    queueSize = 0;
  }
  DBQuery(MYSQL *connp) {
    conn = connp;
    queryStr = "";
    str(queryStr);
  }
  ssize_t execute();
  MYSQL_RES* executeNF();
  std::string preview();
};

void HandleError(MYSQL *conn, DBQuery *query=NULL);


ssize_t LookupReferencePosition(DBQuery &query,
			    std::string &tableName,
			    ssize_t qryPos,
			    ssize_t &refPos);

void GetOneRow(MYSQL *conn, MYSQL_RES *&res, MYSQL_ROW &row);

ssize_t GetNumBlocks(DBQuery &query,
		 std::string &tableName,
		 ssize_t &numBlocks); 

ssize_t GetFirstBlockId(DBQuery &query,
		    std::string &tableName,
		    ssize_t &firstBlock);

ssize_t  GetBlock(DBQuery &query, std::string &tableName,
	      ssize_t refPos, ssize_t qryPos, LAVRow &block);

ssize_t  GetRefBlock(DBQuery &query, std::string &tableName,
		 ssize_t refPos, LAVRow &block);

ssize_t  GetQryBlock(DBQuery &query, std::string &tableName,
		 ssize_t qryPos, LAVRow &block);

ssize_t GetBlock(DBQuery &query,
	     LAVRow &block);

ssize_t GetBlock(MYSQL_ROW &row,
	     LAVRow &block);

ssize_t GetBlock(DBQuery &query,
	     std::string &tableName,
	     ssize_t blockId, 
	     LAVRow &block);

void CreateTable(DBQuery &query,
		 std::string &table,
		 std::vector<std::string> &columns,  
		 std::vector<std::string> &indices,
		 ssize_t dropTable=1);

ssize_t ConvertStr(char *str, ssize_t& value);

ssize_t ConvertStr(char *str, double& value);

ssize_t GetOneFloat(DBQuery &query, double &value);
ssize_t GetOneInt(DBQuery &query, ssize_t &value);
ssize_t GetUpToInt(DBQuery &query, ssize_t &value);

std::string dbstr(std::string &str);

ssize_t FetchRow(MYSQL_ROW &row,
	     LAVRow &block);

ssize_t FetchRow(MYSQL_ROW &row,
	     ssize_t &tStart, ssize_t &tEnd, 
	     ssize_t &qStart, ssize_t &qEnd, ssize_t &strand, ssize_t &chainId);

ssize_t FetchRow(MYSQL_RES *res,
	     ssize_t &tStart, ssize_t &tEnd, 
	     ssize_t &qStart, ssize_t &qEnd, ssize_t &strand, ssize_t &chainId);


void StoreRow(DBQuery &query,
	      std::string &table,
	      ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, 
	      ssize_t strand, ssize_t chainId,
	      ssize_t refSeqId, ssize_t qrySeqId);

void EnqueueStoreRow(DBQuery &query,
		     std::string &table,
		     ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, 
		     ssize_t strand, ssize_t chainId,
		     ssize_t refSeqId, ssize_t qrySeqId, ssize_t maxQueue=-1);

ssize_t ConnectToDB(MYSQL *&conn,
		std::string dbName,
		std::string defaultsFile);

ssize_t ConnectToDB(MYSQL *&conn,
		std::string dbName);

#endif
