/***************************************************************************
 * Title:          NetDB.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _NETDB_H_
#define _NETDB_H_

#include <string>

#include "mysql/mysql.h"
#include "blocks/BlockDB.h"
#include "mctypes.h"

void GetFirstLevelNets(DBQuery &query, std::string &netTableName, 
		       std::vector<ssize_t> &chainIds, ssize_t level = 1, ssize_t orientation = 0);

void BuildChainIdPredicate(std::string &chainTableName, 
			   std::vector<ssize_t> &chainIds,
			   std::string &predicate);

ssize_t CreateTemporaryTable(DBQuery &query,
			 std::string &lavTableName,
			 std::string &chainTableName,
			 std::string &netTableName, 
			 std::string &tempTableName,
			 ssize_t strand,
			 StringSet &tempTableNames,
			 ssize_t maxNumTables, ssize_t verbosity);

ssize_t LookupAlignedBlock(DBQuery &query,
		       std::string &tableName,
		       ssize_t orientation, 
		       ssize_t pos,
		       LAVRow &block,
		       std::set<std::string> &tmpTables, 
		       ssize_t maxNumTables, ssize_t verbosity);

ssize_t LookupAlignedPosition(DBQuery & query,
			  std::string &tableName,
			  ssize_t orientation,
			  ssize_t pos,
			  ssize_t &orthoPos, 
			  std::set<std::string> &tmpTables, ssize_t maxNumTables,
			  ssize_t verbosity);

//
// Lookup the orthologous position of one species in another.
// tableName refers to a table that has a reference seqence aligned to 
// a query sequence.  Pos is a pos in the query sequence.
// orthPos is a position in the reference sequence corresponding to the 
//  query sequence.
// So if there is a position in dog ENm001 that needs to be found in
// human, the table that would be used is human_dog_ENm001. 
//

ssize_t LookupOrthologousPosition(DBQuery & query,
			      std::string &tableName,
			      ssize_t orientation,
			      ssize_t pos,
			      ssize_t refLen,
			      ssize_t &orthPos,
			      ssize_t &orthStrand,
			      std::set<std::string> &tmpTables, ssize_t maxNumTables,
			      ssize_t &deletion,
			      ssize_t skipRevLookup,
			      ssize_t verbosity);


ssize_t FetchSurroundingBlocks(DBQuery &query,
			   std::string &tableName,
			   ssize_t refStartPos, ssize_t refEndPos, 
			   ssize_t qryStartPos, ssize_t qryEndPos,
			   LAVRow &prevRow,
			   LAVRow &nextRow, ssize_t verbosity);


ssize_t GetSequenceLength(DBQuery &query,
		      std::string seqTableName,
		      std::string species,
		      std::string sequence,
		      ssize_t &length);


#endif
