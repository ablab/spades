/***************************************************************************
 * Title:          BlockDB.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <vector>
#include <iostream>

#include "BlockLib.h"
#include "BlockDB.h"
#include "lav/LAVTable.h"


MYSQL_RES* DBQuery::executeNF() {
  // no fault version of execute
  execute();
  MYSQL_RES *res;
  res = mysql_store_result(conn);
  if (res == NULL) {
    std::cout << "error with result: " << mysql_error(conn) << std::endl;
    std::cout << "tried to execute: " << prevQuery << std::endl;
    exit(0);
  }
  return res;
}

ssize_t DBQuery::execute() {
  assert(conn != NULL);
  ssize_t result = mysql_real_query(conn, str().c_str(), str().size());
  if (result) {
    std::cout << "error executing query: " << mysql_error(conn) << std::endl;
    std::cout << "tried to execute: " << str() << std::endl;
  }
  prevQuery = str();
  queryStr = "";
  str(queryStr);
  queueSize = 0;
  return result;
}

std::string DBQuery::preview() {
  return str();
}

void ErrorHandle(MYSQL *conn) {
  std::cout << "mysql error: " << mysql_error(conn) << " " << mysql_errno(conn) << std::endl;
  exit(1);
}


ssize_t LookupReferencePosition(DBQuery &query,
			    std::string &tableName,
			    ssize_t qryPos,
			    ssize_t &refPos) {
  query << "SELECT * FROM " << tableName << " WHERE qStart <= " 
	<< qryPos << " AND qEnd >= " << qryPos;

  query.execute();

  MYSQL_RES *resSet;
  if ((resSet = mysql_store_result(query.conn)) == NULL) {
    std::cout << "getting the orthologous position falied " << std::endl;
    HandleError(query.conn);
    exit(1);
  }

  ssize_t numRows;

  numRows = mysql_num_rows(resSet);
  ssize_t tStart, tEnd, qStart, qEnd;
  refPos = -1;
  if (numRows > 0) {
    MYSQL_ROW row;
    if ((row = mysql_fetch_row(resSet)) == NULL) {
      tStart = atoi(row[1]);
      tEnd   = atoi(row[2]);
      qStart = atoi(row[3]);
      qEnd   = atoi(row[4]);   
      std::cout << "got pos: " << tStart << " " << tEnd << " " << qStart << " " 
		<< qEnd << std::endl;
      GetOffsetPosition(qStart, qEnd, qryPos,
			tStart, tEnd, refPos);

    }		  
  }
}


ssize_t FetchRow(MYSQL_ROW &row,
	     LAVRow &block) {
  return FetchRow(row, block.tStart, block.tEnd, block.qStart, block.qEnd,
		  block.strand, block.chainId);
}

ssize_t  GetBlock(DBQuery &query, std::string &tableName,
	      ssize_t refPos, ssize_t qryPos, LAVRow &block) {
  query << "SELECT * from " << tableName << " where tStart < " << refPos 
	<< " AND tEnd > " << refPos << " AND qStart < " << qryPos 
	<< " AND qEnd > " << qryPos;
  query.execute();
  MYSQL_RES *res = mysql_store_result(query.conn);
  MYSQL_ROW row;
  if (res != NULL and mysql_num_rows(res) > 0) {
    row = mysql_fetch_row(res);
    GetBlock(row, block);
    return 1;
  }
  else {
    return 0;
  }
}




ssize_t  GetRefBlock(DBQuery &query, std::string &tableName,
		 ssize_t refPos, LAVRow &block) {
  query << "SELECT * from " << tableName << " where tStart < " << refPos 
	<< " AND tEnd > " << refPos;
  query.execute();
  MYSQL_RES *res = mysql_store_result(query.conn);
  MYSQL_ROW row;
  if (res != NULL and mysql_num_rows(res) > 0) {
    row = mysql_fetch_row(res);
    GetBlock(row, block);
    return 1;
  }
  else {
    return 0;
  }
}

ssize_t  GetQryBlock(DBQuery &query, std::string &tableName,
		 ssize_t qryPos, LAVRow &block) {
  query << "SELECT * from " << tableName << " where qStart < " << qryPos
	<< " AND qEnd > " << qryPos;
  query.execute();
  MYSQL_RES *res = mysql_store_result(query.conn);
  MYSQL_ROW row;
  if (res != NULL and mysql_num_rows(res) > 0) {
    row = mysql_fetch_row(res);
    GetBlock(row, block);
    return 1;
  }
  else {
    return 0;
  }
}


ssize_t GetBlock(DBQuery &query,
	     LAVRow &block) {
  MYSQL_RES *res;
  MYSQL_ROW row;
  query.execute();
  res = mysql_store_result(query.conn);
  if (res == NULL) 
    return 0;

  row = mysql_fetch_row(res);
  ssize_t retVal;
  if (row != NULL)
    retVal = GetBlock(row, block);
  else 
    retVal = 0;

  mysql_free_result(res);
  return retVal;
}

ssize_t GetBlock(MYSQL_ROW &row,
	     LAVRow &block) {
  return FetchRow(row, block.tStart, block.tEnd, block.qStart, block.qEnd,
		  block.strand, block.chainId);
}

ssize_t GetBlock(DBQuery &query, 
	     std::string &tableName,
	     ssize_t blockId,
	     LAVRow &block) {
  
  query << "SELECT * FROM " << tableName << " WHERE id="<< blockId;
  query.execute();
  MYSQL_RES *res;
  if ((res = mysql_store_result(query.conn))== NULL) HandleError(query.conn);
  
  return FetchRow(res, block.tStart, block.tEnd, block.qStart, block.qEnd,
		  block.strand, block.chainId);

}

void HandleError(MYSQL *conn, DBQuery *query) {
  std::cout << "Mysql error: " << mysql_error(conn) << std::endl;
  if (query != NULL) {
    std::cout << "tried to execute " << query->str() << std::endl;
  }
  exit(0);
}

void GetOneRow(MYSQL *conn, MYSQL_RES *&res, MYSQL_ROW &row) {
  if ((res = mysql_store_result(conn)) == NULL) HandleError(conn);
  if ((row = mysql_fetch_row(res)) == NULL) HandleError(conn);
}
  
ssize_t GetFirstBlockId(DBQuery &query,
		    std::string &tableName,
		    ssize_t &firstBlock) {
  query << "select min(id) from " << tableName;
  query.execute();
  MYSQL_RES *res;
  MYSQL_ROW row;
  GetOneRow(query.conn, res, row);
  firstBlock = atoi(row[0]);
}

std::string dbstr(std::string &str)
{
  return "\"" + str + "\"";
}

ssize_t GetNumBlocks(DBQuery &query,
		 std::string &tableName, 
		 ssize_t &numBlocks) { 
  query << "select max(id) from " << tableName;
  query.execute();
  MYSQL_RES *res;
  MYSQL_ROW row;
  GetOneRow(query.conn, res, row);
  numBlocks = atoi(row[0]);
  mysql_free_result(res);
}

ssize_t GetUpToInt(DBQuery &query, ssize_t &value) {
 query.execute(); 

 MYSQL_RES *res; 
 if ((res = mysql_store_result(query.conn)) == NULL) HandleError(query.conn);
 ssize_t numRows;
 numRows = mysql_num_rows(res);
 if (numRows <= 0) {
   return 0;
 }
 MYSQL_ROW row; 
 if (numRows > 0) {
   if ((row = mysql_fetch_row(res)) == NULL) {
     std::cout <<"error: " << mysql_error(query.conn) << std::endl; 
     exit(1); 
   } 
   if (row[0] != NULL)
     ConvertStr(row[0], value); 
 }
 mysql_free_result(res);
 return 1;
}

ssize_t GetOneFloat(DBQuery &query, double &value) {
  // Get an int, and say there is an error with query.
  
  query.execute(); 

  MYSQL_RES *res; 
  if ((res = mysql_store_result(query.conn)) == NULL) HandleError(query.conn);
  if (mysql_num_rows(res) == 0) { 
    std::cout << "there must be at least an ssize_t " << std::endl; 
    std::cout << "failed query: " << query.prevQuery << std::endl;
    exit(0);
  }
  MYSQL_ROW row; 
  if ((row = mysql_fetch_row(res)) == NULL) {
    std::cout <<"error: " << mysql_error(query.conn) << std::endl; 
    exit(1); 
  } 
  if (row[0] != NULL)
    ConvertStr(row[0], value); 
  mysql_free_result(res);
}


ssize_t GetOneInt(DBQuery &query, ssize_t &value) {
  // Get an int, and say there is an error with query.
  
  query.execute(); 

  MYSQL_RES *res; 
  if ((res = mysql_store_result(query.conn)) == NULL) HandleError(query.conn);
  if (mysql_num_rows(res) == 0) { 
    std::cout << "there must be at least an ssize_t " << std::endl; 
    std::cout << "failed query: " << query.prevQuery << std::endl;
    exit(0);
  }
  MYSQL_ROW row; 
  if ((row = mysql_fetch_row(res)) == NULL) {
    std::cout <<"error: " << mysql_error(query.conn) << std::endl; 
    exit(1); 
  } 
  if (row[0] != NULL)
    ConvertStr(row[0], value); 
  mysql_free_result(res);
  return 1;
}

inline ssize_t ConvertStr(char *str, ssize_t& value) { value = atoi(str); return 1; }

inline ssize_t ConvertStr(char *str, double& value) { value = atof(str); return 1; }


void CreateTable(DBQuery &query,
		 std::string &table,
		 std::vector<std::string> &columns,  
		 std::vector<std::string> &indices,
		 ssize_t drop) {
  if (drop){
    query << "DROP TABLE IF EXISTS " << table << "; ";
    query.execute();
  }
  query << "CREATE TABLE IF NOT EXISTS " << table ;

  ssize_t i;
  // sanity check
  assert(columns.size() >= 2 || 
	 (printf("columnsize must be greater than 2") == 0));
  
  query << "(";
  for (i = 0; i < columns.size() - 2; i+= 2) 
    query << " " << columns[i] << " " << columns[i+1] << ", ";
  query << " " << columns[columns.size()-2] << " " 
	       << columns[columns.size()-1];
  
  if (indices.size() > 0) {
    query << ", ";
    for (i = 0; i < indices.size()-1; i++) 
      query << " INDEX ( " << indices[i] << "), ";
    query << " INDEX ( " << indices[i] << ") ";
  }
    
  query << ")";
//  std::cout << "creating table with: " << query.str() << std::endl;
  if (query.execute()) 
    HandleError(query.conn);
}


void EnqueueStoreRow(DBQuery &query,
		     std::string &table,
		     ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, 
		     ssize_t strand, ssize_t chainId,
		     ssize_t refSeqId, ssize_t qrySeqId, ssize_t maxQueue) {
  query.queueSize++;
  if (query.str() == "") {
  query << "INSERT INTO " << table 
	<< " (tStart, tEnd, qStart, qEnd, strand, chainId, tid, qid) "
	<< " VALUES ( " 
	<< tStart << "," 
	<< tEnd << "," 
	<< qStart << "," 	
	<< qEnd << "," 
	<< strand <<  ","
	<< chainId << ","
	<< refSeqId << "," 
	<< qrySeqId << ")";
  }
  else {
    query << ",( " 
	  << tStart << "," 
	  << tEnd << "," 
	  << qStart << "," 	
	  << qEnd << "," 
	  << strand <<  ","
	  << chainId << ","
	  << refSeqId << "," 
	  << qrySeqId << ")";
  }
  if (maxQueue > 0 && query.queueSize >= maxQueue) {
    if (query.execute())
      HandleError(query.conn);
  }
}

void StoreRow(DBQuery &query,
	      std::string &table,
	      ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, 
	      ssize_t strand, ssize_t chainId,
	      ssize_t refSeqId, ssize_t qrySeqId) {

  query << "INSERT INTO " << table 
	<< " (tStart, tEnd, qStart, qEnd, strand, chainId, tid, qid) "
	<< " VALUES ( " 
	<< tStart << "," 
	<< tEnd << "," 
	<< qStart << "," 	
	<< qEnd << "," 
	<< strand <<  ","
	<< chainId << ","
	<< refSeqId << "," 
	<< qrySeqId << ")";

  if (query.execute()) 
    HandleError(query.conn);
    
}

ssize_t ConnectToDB(MYSQL *&conn,
		std::string dbName) {
  ConnectToDB(conn, dbName, "");
}

ssize_t ConnectToDB(MYSQL *&conn,
		std::string dbName,
		std::string defaultsFile) {

  conn = mysql_init(NULL);
  if (defaultsFile != "") {
    mysql_options(conn, MYSQL_READ_DEFAULT_FILE, defaultsFile.c_str());
    mysql_options(conn, MYSQL_READ_DEFAULT_GROUP, "client");
  }
  if (conn == NULL) {
    std::cout << "Could not initialize mysql connection. " << std::endl;
    exit(1);
  }
  //  if ((conn = mysql_real_connect(conn, NULL, "mchaisso", NULL, NULL, 0, "/scratch/mchaisso/mysql/global/mysql.socket", 0)) 
	char *user = getenv("USER");
	if (user == NULL) {
		std::cout << "ERROR, must 'USER' must be set to the user running this" << std::endl;
		std::cout << "program. " << std::endl;
		exit(1);
	}
  if ((conn = mysql_real_connect(conn, NULL, user, NULL, NULL, 0, NULL, 0)) 
      == NULL) {
    if (conn != NULL)
      HandleError(conn);
    else {
      std::cout << "Could not connect to mysql. " << std::endl;
      exit(1);
    }
  }

  if (mysql_select_db(conn, dbName.c_str()))
    HandleError(conn);
  return 1;
}

ssize_t FetchRow(MYSQL_ROW &row,
	     ssize_t &tStart, ssize_t &tEnd, 
	     ssize_t &qStart, ssize_t &qEnd, ssize_t &strand, ssize_t &chainId) {
  tStart = atoi(row[1]);
  tEnd   = atoi(row[2]);
  qStart = atoi(row[3]);
  qEnd   = atoi(row[4]);
  strand = atoi(row[5]);
  chainId = atoi(row[6]);
  
  return 1;
}

ssize_t FetchRow(MYSQL_RES *res,
	     ssize_t &tStart, ssize_t &tEnd, 
	     ssize_t &qStart, ssize_t &qEnd, ssize_t &strand, ssize_t &chainId) {
  MYSQL_ROW row;
  if (! mysql_num_rows(res) > 0) {
    return 0;
  }
  if ((row = mysql_fetch_row(res)) == NULL) {
    std::cout << "erro fetchign row:" << std::endl;
    return 0;
  }
  return FetchRow(row, tStart, tEnd, qStart, qEnd, strand, chainId);
}
