/***************************************************************************
 * Title:          NetDB.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "NetDB.h"
#include "blocks/BlockLib.h"

#include <sstream>

void GetFirstLevelNets(DBQuery &query, std::string &netTableName, 
		       std::vector<ssize_t> &chainIds, ssize_t level, ssize_t dir) {

  // Find all parts of a net that aren't overlapping, and are on the 
  query << "SELECT chainid FROM " << netTableName << " WHERE level = " << level 
	<< " AND orientation =" << dir;
  query.execute();

  MYSQL_RES *res;
  MYSQL_ROW row;

  res = mysql_store_result(query.conn);
  
  while ((row = mysql_fetch_row(res)) != NULL) {
    chainIds.push_back(atoi(row[0]));
  }

  mysql_free_result(res);
}


void BuildChainIdPredicate(std::string &chainTableName, 
			   std::vector<ssize_t> &chainIds,
			   std::string &predicate) {

  ssize_t i;
  std::stringstream predstr;
  if (chainIds.size() > 0) {
    predstr << " (";
    for (i = 0; i < chainIds.size() - 1; i++) {
      predstr << chainTableName << ".chainId = " << chainIds[i] << " OR ";
    }
    predstr << chainTableName << ".chainId = " << chainIds[i] << " ) ";
  }
  predicate = predstr.str();
}

ssize_t CreateTemporaryTable(DBQuery &query,
			 std::string &lavTableName,
			 std::string &chainTableName,
			 std::string &netTableName, 
			 std::string &tempTableName,
			 ssize_t strand,
			 StringSet &tempTableNames,
			 ssize_t maxNumTables, 
			 ssize_t verbosity) {
  // Ask the database if the temporary table already exists
  return 1;
  query << "show tables like \""<< tempTableName << "\"";
  //  query << "select min(tStart) from " << tempTableName;

  if (!query.execute()) {
    MYSQL_RES *res;
    MYSQL_ROW row;
    res = mysql_store_result(query.conn);
    if (mysql_num_rows(res) > 0) {
      //std::cout << tempTableName << " already exists, not rebuilding " << std::endl;
      return 1;  // the table already exists
    }
  }


  // Create a table of lav alignments that are chained together.
  if (tempTableNames.find(tempTableName) == tempTableNames.end()) {
    // Look to see fi there are too many tables
    StringSet::iterator removeIt;
    if (tempTableNames.size() > maxNumTables) {
      // Too many tables. Eject one.
      std::string goodbye;
      removeIt = tempTableNames.begin();
      goodbye = *removeIt;
      query << "DROP TABLE " << goodbye;
      query.execute();
      //      std::cout << "got rid of " << *removeIt << std::endl;
      tempTableNames.erase(removeIt);
    }
	      
	      
    std::vector<ssize_t> firstLevel;
    std::string firstLevelP = "";
    ssize_t level = 1;
    if (strand == 1) {
      level = 3;
    }
    GetFirstLevelNets(query, netTableName, firstLevel, level , strand );
    BuildChainIdPredicate(chainTableName, firstLevel, firstLevelP);
    if (firstLevelP != "") {
      query << "CREATE TEMPORARY TABLE " << tempTableName << " AS SELECT " 
	    << lavTableName << ".id, " <<   lavTableName << ".tStart,"
	    << lavTableName << ".tEnd," <<  lavTableName << ".qStart, "
	    << lavTableName << ".qEnd," <<  lavTableName << ".strand, " 
	    << lavTableName << ".chainId FROM " << lavTableName << ", " 
	    << chainTableName << " WHERE " << lavTableName << ".chainId = " 
	    << chainTableName << ".lavId AND " << lavTableName << ".strand =" << strand
	    << " AND " << firstLevelP;
      if (verbosity > 0 ) {
	std::cout << "creating temporary table: " << tempTableName << " with: ";
	std::cout << query.str() << std::endl;
      }
      query.execute();
      tempTableNames.insert(tempTableName);
      return 1;
    }
  }
  return 0;
}

ssize_t LookupAlignedBlock(DBQuery &query,
		       std::string &tableName,
		       ssize_t orientation, 
		       ssize_t pos,
		       LAVRow &block,
		       std::set<std::string> &tmpTables, 
		       ssize_t maxNumTables, ssize_t verbosity) {
  MYSQL_RES *res;
  MYSQL_ROW row;

  ssize_t maxScore, maxChainId;
  std::vector<ssize_t> chainIds;
  std::string predicate;

  std::string chainTableName;
  std::string netTableName;
  std::string tempTableName;
  chainTableName = tableName + "_chain";
  netTableName   = tableName + "_net";
  tempTableName   = tableName + "_temp";

  CreateTemporaryTable(query, 
                       tableName, chainTableName, netTableName, tempTableName,
		       orientation,
		       tmpTables, maxNumTables, verbosity);
    
  // Find the coordinate in the new/target alignment
  query << "SELECT id, tStart, tEnd, qStart,qEnd,strand,chainID FROM "
	<< tempTableName << " WHERE qStart <= " << pos << " AND qEnd >= " 
	<< pos << " AND strand = " << orientation;  
  if (verbosity > 0 ) {
    std::cout << "getting aligned block with: " << query.str() << std::endl;
  }
  return GetBlock(query, block);
}


ssize_t LookupAlignedPosition(DBQuery & query,
			  std::string &tableName,
			  ssize_t orientation, ssize_t pos,
			  ssize_t &orthPos,
			  std::set<std::string> &tmpTables, 
			  ssize_t maxNumTables, ssize_t verbosity) {
  LAVRow row;
  if (LookupAlignedBlock(query, tableName, orientation, pos, 
			 row,
			 tmpTables, maxNumTables, verbosity)) {
    if (verbosity > 0) {
      std::cout << " found aligned block " << std::endl;
    }
    GetOffsetPosition(row.qStart, row.qEnd, pos +1,
		      row.tStart, row.tEnd, orthPos);
    return 1;
  }
  return 0;
}

ssize_t LookupOrthologousPosition(DBQuery & query,
			      std::string &tableName,
			      ssize_t orientation, ssize_t pos,
			      ssize_t refLen,
			      ssize_t &orthPos,
			      ssize_t &orthStrand,
			      std::set<std::string> &tmpTables, ssize_t maxTempTables,
			      ssize_t &deleted,
			      ssize_t skipRevLookup,
			      ssize_t verbosity) {

  MYSQL_RES *res;
  MYSQL_ROW row;
  ssize_t rowFound;
  // Look for the position in a forward strand alignment on the main alignment
  std::string chainTableName;
  std::string netTableName;
  std::string tempTableName;
  chainTableName = tableName + "_chain";
  netTableName   = tableName + "_net";
  tempTableName   = tableName + "_temp";
  deleted = 0;
  orthPos = -1;
  rowFound = LookupAlignedPosition(query, 
				   tableName, 
				   orientation, pos,
				   orthPos, tmpTables, maxTempTables, verbosity);
  if (rowFound) {
    orthStrand = 0;
    return 1;
  }
  

  // Didn't find it there, look through reverse strand alignments.
  if (! skipRevLookup ) {
    ssize_t revPos;
    revPos = refLen - pos + 1; 
    ssize_t strand = 1;
    if (orientation == 1) 
      strand = 0;
    else
      strand = 1;
    
    LAVRow revStrandRow;
    //  std::cout << "creating temporay table " << std::endl;
    std::string revTempTableName = tempTableName + "_reverse";
    if (CreateTemporaryTable(query, tableName, chainTableName, netTableName, revTempTableName,
			     !orientation,
			     tmpTables, maxTempTables, verbosity)) {
      //  std::cout << "done creating temporary table " << std::endl;
      query << "SELECT * FROM "
	    << revTempTableName << " WHERE qStart <= " << revPos << " and qEnd >= " 
	    << revPos << " AND strand = " << strand; 
      if ( verbosity > 0 ) {
	std::cout << "now looking up " << query.str() << std::endl;
      }
      if (GetBlock(query, revStrandRow)) {
	GetOffsetPosition(revStrandRow.qStart, revStrandRow.qEnd, revPos,
			  revStrandRow.tStart, revStrandRow.tEnd, orthPos);
	orthStrand = 0;
	return 1;
      }
    }
  }
  
  // Didn't find an alignment in the opposite orientation that works.  
  // Look between two alignments and extrapolate the position.
  LAVRow prevRow, nextRow;
  FetchSurroundingBlocks(query, tempTableName,
			 -1, -1,
			 pos, pos,
			 prevRow, nextRow, verbosity);
  
  if (nextRow.qStart == -1 and 
      prevRow.qStart == -1) {
    // Something bad here.  Didn't find an surrounding blocks.  Perhaps return 0.
    return 0;
  }

  ssize_t maxQryStart, maxRefStart, minQryEnd, minRefEnd;
  ssize_t qryStart, refStart, qryEnd, refEnd;
  if (nextRow.qStart != -1) {
    //    maxQryStart = nextRow.qStart;
    //    maxRefStart = nextRow.tStart;
    qryEnd = nextRow.qStart;
    refEnd = nextRow.tStart;
  }
  else {
    // No ending interval was found.  This means that the position maps 
    // past the end of the human sequence.
    //maxQryStart = pos;
    qryEnd = pos;
    refEnd = prevRow.tEnd + (pos - prevRow.qEnd);
  } 

  if (prevRow.qStart != -1) {
    //    minQryEnd = prevRow.qEnd;
    //    minRefEnd = prevRow.tEnd;
    qryStart = prevRow.qEnd;
    refStart = prevRow.tEnd;
  }
  else {
    // No beginning interval was found.
    // The sequence maps before the first alignment in human.
    qryStart = pos; //minQryEnd   = pos;
    refStart = nextRow.tStart  - (nextRow.qStart - pos);
  }

  // now extrapolate between rows.
  ssize_t qryGapLen;

  if ( ((double( refEnd - refStart ) / ( qryEnd - qryStart ))  > 10) or 
       ((double( qryEnd - qryStart ) / ( refEnd - refStart ))  > 10))
    deleted = 1;
      
  GetOffsetPosition(qryStart, qryEnd, pos,
		    refStart, refEnd, orthPos);


  // Since this is in an unaligned region, it's possible that the
  // sequence is simply not known. 

  return 1;
}
  
ssize_t FetchSurroundingBlocks(DBQuery &query,
			   std::string &tableName,
			   ssize_t refStartPos, ssize_t refEndPos, 
			   ssize_t qryStartPos, ssize_t qryEndPos,
			   LAVRow &prevRow,
			   LAVRow &nextRow, ssize_t verbosity) {

  std::string endCol, startCol;
  ssize_t startPos, endPos;
  if (refStartPos < 0) {
    endCol   = "qEnd";
    startCol = "qStart";
    startPos = qryStartPos;
    endPos   = qryEndPos;
  }
  else {
    endCol   = "tEnd";
    startCol = "tStart";
    startPos   = refStartPos;
    endPos     = refEndPos;
  }
  // Still didn't find an alignment containing this position.

  // Find the endpoints of the alignments surrounding this position.  
  // The left endpoint is maxqend, right minqstart.

  query << "SELECT max(" << endCol << ") FROM " << tableName 
	<< " WHERE " << startCol << " <= " << startPos;
  if (verbosity > 1) {
    std::cout << "getting surrounding block with: " << query.str() << std::endl;
  }
  ssize_t maxEnd, minStart; 

  maxEnd = -1; minStart = -1;

  if (GetUpToInt(query, maxEnd) == 0) {
    std::cout << "didn't get anything from: " << query.prevQuery << std::endl;
    return 0;
  }

  query << "SELECT min(" << startCol << ") FROM " << tableName 
	<< " WHERE " << endCol << " >= " << endPos;
  if (verbosity > 1) {
    std::cout << "getting end col with " << query.str() << std::endl;
  }
  if  (GetUpToInt(query, minStart) == 0) {
    std::cout << "didn't get anything from: " << query.prevQuery << std::endl;
    return 0;
  }

  if (minStart != -1) {
    query << "SELECT * FROM " << tableName 
	  << " WHERE " << startCol << " = " << minStart;
    if  (!GetBlock(query, prevRow))
      return 0;

    // Now check to see if this row includes this point.
    ssize_t diff;
    if (refStartPos >= 0) {
      if (prevRow.tEnd > refStartPos) {
	diff = prevRow.tEnd - refStartPos;
	prevRow.tEnd = refStartPos -1;
	prevRow.qEnd = prevRow.qEnd - diff - 1;
      }
    }
    else {
      if (prevRow.qEnd > qryStartPos) {
	diff = prevRow.qEnd - qryStartPos;
	prevRow.tEnd = prevRow.tEnd - diff - 1;
	prevRow.qEnd = qryStartPos - 1;
      }
    }
  }
  ssize_t diff;
  if (maxEnd != -1) {
    query << "SELECT * FROM " 
	  << tableName << " WHERE " << endCol << " = " << maxEnd ;
  
    if (!GetBlock(query, nextRow))
      return 0;
    if (refEndPos >= 0) {
      if (nextRow.tStart < refEndPos) {
	diff = nextRow.tStart - refEndPos;
	nextRow.tStart = refEndPos + 1;
	nextRow.qStart = nextRow.qStart + diff + 1;
      }
    }
    else {
      if (nextRow.qStart < qryEndPos) {
	diff = nextRow.qStart - qryEndPos;
	nextRow.tStart = nextRow.tStart + diff;
	nextRow.qStart = qryEndPos + 1;
      }
    }

  }
  return 1;
}

ssize_t GetSequenceLength(DBQuery &query,
		      std::string seqTableName,
		      std::string species,
		      std::string sequence,
		      ssize_t &length) {
  query << "SELECT length FROM " << seqTableName << " WHERE name =" << dbstr(species) 
	<< " AND sequence="<< dbstr(sequence) ;
  if (GetOneInt(query, length))
    return 1;
  else 
    return 0;
}

