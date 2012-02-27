/***************************************************************************
 * Title:          ChainToDB.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <unistd.h>
#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include "mysql/mysql.h"

#include "chain/ChainReader.h"
#include "chain/Chain.h"
#include "blocks/BlockDB.h"
#include "blocks/BlockLib.h"
#include "lav/LAVTable.h"
#include "utils.h"

void StoreIds(DBQuery &query, std::string &table, 
	      std::vector<ssize_t> &lavIds, std::vector<ssize_t> &chainIds);

void EnqueueStore(DBQuery &query, std::string &table, 
		  std::vector<ssize_t> &lavIds, std::vector<ssize_t> &chainIds, 
		  ssize_t lavId, ssize_t chainId, ssize_t maxQueue);

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &defaultsFile,
	    std::string &refName,
	    std::string &qryName,
	    std::string &chaindir,
	    std::string &lavedir,
	    std::string &seqName,
	    std::string &database,
	    std::string &chainFileName,
	    ssize_t &dropTable,
	    ssize_t &printTabular);


void PrintUsage();


int main(int argc, char* argv[]) {
  std::string  refName;
  std::string  qryName;
  std::string  chainDir;
  std::string  lavDir;
  std::string  seqName;
  std::string  lavTable;
  std::string chainTable;
  std::string database;
  std::string tabularName;
  std::string defaultsFile;
  ssize_t dropTable = 0;
  ssize_t printTabular = 0;

  chainDir = "bzchain";
  lavDir   = "bzalign";
  tabularName = "";
  std::cout << "starting chaintodb " << std::endl;
  std::string chainFileName;

  chainFileName = "";
  InitEnv(argc, argv, defaultsFile, 
	  refName, qryName, seqName, 
	  chainDir, lavDir, database, chainFileName, dropTable, printTabular);
  
  std::string lavDBName;
  

  std::vector<Chain*> chains;
  
  std::ofstream tabOutFile;
  // Read the chain file
  if (chainFileName == "") {
    chainFileName = chainDir + "/" + refName +  "."
      + qryName + ".chain";
  }
  tabularName = chainFileName + ".txt";

  if (printTabular) {
    openck(tabularName, tabOutFile, std::ios::out);
  }

  ChainReader::ReadChainFile(chainFileName, chains);
  std::cout << "read the chain file " << std::endl;
  // Initialize the connection to the alignments

  
  // Create the lav table name
  lavTable = refName + "_" + qryName + "_" + seqName;

  // Create the chain table 
  chainTable = refName + "_" + qryName + "_" + seqName + "_chain";

  MYSQL *conn;
  ConnectToDB(conn, database.c_str(), defaultsFile);

  DBQuery query;
  query.conn = conn;
  
  if (dropTable) {
    query << "DROP TABLE IF EXISTS " << chainTable;
    query.execute();
  }

  // Chain table is a series of references to alignments stored in
  // the lav file.

  std::vector<std::string> columns, indices;
  
  columns.push_back("id"); columns.push_back(" ssize_t AUTO_INCREMENT PRIMARY KEY ");
  columns.push_back("chainid"); columns.push_back(" ssize_t ");
  columns.push_back("lavid"); columns.push_back(" ssize_t ");
  
  indices.push_back("chainid");
  indices.push_back("lavid");
  
  std::cout << "creating " << chainTable << std::endl;
  CreateTable(query, chainTable, columns, indices);
  

  // The chain format stores each alignment as a series of ungapped
  // alignments.  Each ungapped alignment *should* correspond to
  // one of the blastz alignments.  I'm not sure if the chaining
  // fills in extra gaps or not.  I don't think it does, but I'm not sure.
  // I guess now's the time to find out.
  MYSQL_RES *res;
  MYSQL_ROW row;

  ssize_t c, a; // iterate over chains, alignments
  Chain *chain;
  ssize_t tPos, qPos;
  ssize_t tEnd, qEnd;
  LAVRow lavRow;
  std::vector<ssize_t> chainIds;
  std::vector<ssize_t> lavIds;
  std::string tempLavTable;
  tempLavTable = lavTable + "_tmp";
  std::cout << "creating temporary lav table: " << tempLavTable << std::endl;
  query << "CREATE TEMPORARY TABLE " << tempLavTable << " AS select * from " << lavTable;
  query.execute();
  query << "ALTER TABLE " << tempLavTable << " TYPE=HEAP ";
  query.execute();
  query << "ALTER TABLE " << tempLavTable << " ADD INDEX (tStart) ";
  query.execute();
  query << "ALTER TABLE " << tempLavTable << " ADD INDEX (tEnd) ";
  query.execute();
  query << "ALTER TABLE " << tempLavTable << " ADD INDEX (qStart) ";
  query.execute();
  query << "ALTER TABLE " << tempLavTable << " ADD INDEX (qEnd) ";
  query.execute();


  //  query << "LOCK TABLES " << chainTable << " READ " ;
  //  std::cout << query.str() << std::endl;
  //  query.execute();
  ssize_t tot = 0;
  for (c = 0; c < chains.size(); c++) 
    tot += chains[c]->numAlign();
  std::cout << "processing : " << tot << " chains " << std::endl;
  ssize_t notFoundRows = 0;
  ssize_t tooManyRows = 0;
  ssize_t storedRows = 0;
      
  ssize_t t = 0;
  for (c = 0; c < chains.size(); c++) {
    chain = chains[c];
    tPos = chain->header.tStart;
    qPos = chain->header.qStart;
    for (a = 0; a < chain->numAlign(); a++) {
      if (t % 1000 == 0) 
	std::cout << t << " / " << tot << std::endl;
      ++t;
      tEnd = tPos + chain->size[a];
      qEnd = qPos + chain->size[a];
      
      // look up the coordinates in the database
      /*
	std::cout << "looking up " << tempLavTable << " " << tPos + 1 << " " << tEnd
	<< " " << qPos + 1 << " " << qEnd << std::endl;
      */
      query << "SELECT * from " << tempLavTable << " WHERE " 
	    << " tStart = " << tPos + 1 << " AND tEnd = " << tEnd 
	    << " AND qStart = " << qPos + 1 << " AND qEnd = " << qEnd;

      if (query.execute()) {
	HandleError(query.conn);
      }

      res = mysql_store_result(query.conn);
      
      if (res != NULL && mysql_num_rows(res) == 0) {
	mysql_free_result(res);
	query << "SELECT * from " << tempLavTable << " WHERE " 
	      << " tStart <= " << tPos + 1 << " AND tEnd >= " << tEnd 
	      << " AND qStart <= " << qPos + 1 << " AND qEnd >= " << qEnd << 
	  " order by tStart - tEnd";
	if (query.execute()) {
	  HandleError(query.conn);
	}
	res = mysql_store_result(query.conn);
      }
      if (res != NULL && mysql_num_rows(res) > 0) {
	row = mysql_fetch_row(res);
	GetBlock(row, lavRow);
	if (mysql_num_rows(res) > 1) {
	  ++tooManyRows;
	  /*	  std::cout << " chain: " << c << " , " << a 
		  << " has too many rows : " << lavRow.chainId 
		  << " : " << mysql_num_rows(res) << std::endl;
		  std::cout << "coords: " << tPos +1 << " " << tEnd 
		  << " " << qPos + 1 << " " <<  qEnd << std::endl;
	  */
	}
	else { 
	  if (!printTabular) {
	    EnqueueStore(query, chainTable, lavIds, chainIds, lavRow.chainId, chain->id, 100);
	  }
	  else {
	    tabOutFile << 0 << "\t" << lavRow.chainId << "\t" << chain->id << std::endl;
	  }
	  ++storedRows;
	  //	  std::cout << "found it " << std::endl;
	}
	mysql_free_result(res);
      }
      else {
	++notFoundRows;
	/*
	  std::cout << query.prevQuery << " did not return any rows at " 
		  << c << ", " << a << std::endl;
	*/
      }
      tPos += (chain->size[a] + chain->dt[a]);
      qPos += (chain->size[a] + chain->dq[a]);
      
    } 
  }  
  std::cout << "num chains " << chains.size() << " stored: " << storedRows 
	    << " multiple row matched: " << tooManyRows 
	    << " rows not found: " << notFoundRows << std::endl;
  if (!printTabular) 
    StoreIds(query, chainTable, lavIds, chainIds);
  else
    tabOutFile.close();
  //  query << "UNLOCK TABLES ";
  //  query.execute();
}

void StoreIds(DBQuery &query, std::string &table, 
	      std::vector<ssize_t> &lavIds, std::vector<ssize_t> &chainIds) {
  if (lavIds.size() == 0) 
    return;

  query << " INSERT INTO " << table << "(chainid, lavid) VALUES (" << chainIds[0]
	<< "," << lavIds[0] << ")";
  ssize_t i;
  for (i = 1; i < lavIds.size(); i++) {
    query << ", ( " << chainIds[i] << ", " << lavIds[i] << ")";
  }
  //  std::cout << "executing " << query.str() << std::endl;
  query.execute();
  lavIds.clear();
  chainIds.clear();
}


void EnqueueStore(DBQuery &query, std::string &table, 
		  std::vector<ssize_t> &lavIds, std::vector<ssize_t> &chainIds, 
		  ssize_t lavId, ssize_t chainId, ssize_t maxQueue) {
  lavIds.push_back(lavId);
  chainIds.push_back(chainId);

  if (lavIds.size() >= maxQueue) {
    StoreIds(query, table, lavIds, chainIds);
  }
}

ssize_t InitEnv(int argc, char* argv[], 
	    std::string &defaultsFile,
	    std::string &refName,
	    std::string &qryName,
	    std::string &seqName,
	    std::string &chainDir,
	    std::string &lavDir,
	    std::string &database,
	    std::string &chainFileName,
	    ssize_t &dropTable,
	    ssize_t &printTabular) {
  defaultsFile = "";

  ssize_t copt;
  // no options for now
  while ( (copt=getopt(argc, argv, "f:c:l:dtC:")) != EOF){
    switch(copt) {
    case 'f':
      defaultsFile = optarg;
      continue;
    case 'd':
      dropTable = 1;
      continue;
    case 'C':
      chainFileName = optarg;
      continue;
    case 'c':
      chainDir = optarg;
      continue;
    case 'l':
      lavDir = optarg;
      continue;
    case 't':
      printTabular = 1;
      continue;
    default:
      PrintUsage();
      exit(0);
    }
  }
  ssize_t ind = optind;

  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  database = argv[ind];
  ind++;

  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  refName = argv[ind];
  ind++;
  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  qryName = argv[ind];
  ind++;
  if (ind >= argc) {
    PrintUsage();
    exit(0);
  }
  seqName = argv[ind];
}


void PrintUsage() 
{
  std::cout  << "usage: chainToDB [-c chain_dir (bzchain)] [-f defaults-file] "
	     << " [-C chainFile ] "  
	     << " database refSpec qrySpec seqName  "
	     << "[-l lavdir (bzalign)] " << std::endl;
}

