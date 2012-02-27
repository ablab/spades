/***************************************************************************
 * Title:          FindConnections.cpp 
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
#include "dbOrthoPos.h"
#include "InversionBins.h"
#include "InversionUtils.h"


int main(int argc, char*argv[]) {

  std::string dbName, seqName, binFileName, outFileName;
    
  if (argc <= 4) {
    std::cout << "usage: findconnections dbname seqname binfile outfile " << std::endl;
    exit(1);
  }
  
  dbName      = argv[1];
  seqName     = argv[2];
  binFileName = argv[3];
  outFileName = argv[4];
  
  StringVector species;
  std::vector<ssize_t> humanStartPos,  humanEndPos;
  BinMap binnedInversions;
  
  // Connect to the database;
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);
  // Create the query and associate it with the databae
  DBQuery query(dbConn);

  ReadBinFile(binFileName, species, binnedInversions);
  

  BinMap::iterator binIt;
  StringVector::iterator specIt, querySpecIt;
  ssize_t specIndex;
  ssize_t specStartPos, specEndPos;
  StringSet tempTableNames;
  std::string tempTableName, netTableName, chainTableName, lavTableName;
  ssize_t maxTempTables;


  ssize_t seqLength, seqStartPos, seqEndPos, seqStartStrand, seqEndStrand;
  ssize_t seqStartDeleted, seqEndDeleted;

  maxTempTables = 1024;
  std::map<std::string, ssize_t> specStart, specEnd;
  std::map<std::string, DNASequence *> seqMap;
  ssize_t binNumber = 0;
  for (binIt = binnedInversions.begin();
       binIt != binnedInversions.end();
       ++binIt) {
    // Each bin represents a separate ortholgous inversion. 
    // There is one graph per inversion.
    // Look through every species to see if it is listed in the bin or not.
    // If it is, we can just use the values in the bin for constructing the graph
    // It it isn't, add another entry to the bin.

    specStart.clear();
    specEnd.clear();

    // Locate the starting and ending position of the inversion for each species.
    for (specIt = species.begin(); specIt != species.end(); ++specIt) {
      specIndex = FindSpecies((*binIt).second, *specIt);
      if (specIndex != -1) {
	specStart[*specIt] = (*binIt).second[specIndex]->startPos;
	specEnd[*specIt] = (*binIt).second[specIndex]->endPos;
      }
      else {
	lavTableName  = *specIt + "_human_" + seqName;
	tempTableName = lavTableName + "_temp";
	netTableName  = lavTableName + "_net";
	chainTableName= lavTableName + "_chain";
	
	CreateTemporaryTable(query, lavTableName, chainTableName, netTableName, tempTableName, 
			     tempTableNames, maxTempTables);
	GetSequenceLength(query, "sequences", *specIt, seqName, seqLength);
	// (*binIt).startPos is the position in human.
	// seqStartPos is the starting position of the inversion in the other sequence.
	std::cout << "looking up pos for " << *specIt << std::endl;
	seqStartPos = -1;
	seqEndPos   = -1;
	LookupOrthologousPosition(query, lavTableName, 0, humanStartPos[binNumber], seqLength, 
				  seqStartPos, seqStartStrand, tempTableNames, maxTempTables, seqStartDeleted);

	LookupOrthologousPosition(query, lavTableName, seqStartStrand, humanEndPos[binNumber], seqLength, 
				  seqEndPos, seqEndStrand, tempTableNames, maxTempTables, seqEndDeleted);

	if (seqStartDeleted and seqEndDeleted) {
	  specStart[*specIt] = -1;
	  specEnd[*specIt] = -1;
	}
	else {
	  if (seqStartPos < seqEndPos) {
	    specStart[*specIt] = seqStartPos;
	    specEnd[*specIt] = seqEndPos;
	  }
	  else {
	    specStart[*specIt] = seqEndPos;
	    specEnd[*specIt] = seqStartPos;
	  }
	  std::cout << " for species: " << *specIt << " got: " << seqStartPos 
		    << " " << seqEndPos << std::endl;
	}
	std::cout << "found human pos: " << *specIt<< " " 
		  << humanStartPos[binNumber] << " " << humanEndPos[binNumber] << " "
		  << seqStartPos << " " << seqEndPos << std::endl;
      }
    }
    ssize_t validation;
    for (specIt = species.begin(); specIt != species.end(); ++specIt) {
      // Look to see if there is already an entry at this inversion for this species.
      specIndex = FindSpecies((*binIt).second, *specIt);
      
      if (specIndex == -1) {
	// There is no entry for this.  Validate the inversion of this sequence in all others.
	std::cout << "No entry found for " << *specIt << " " 
		  << humanStartPos[binNumber] << " " << humanEndPos[binNumber] << " checking:";
	
	ValidatedInversion *valInv;
	if (specStart[*specIt] != -1 and specEnd[*specIt] != -1) {
	  valInv = new ValidatedInversion;
	  valInv->startPos = specStart[*specIt];
	  valInv->endPos   = specEnd[*specIt];
	  valInv->species  = *specIt;
	  for (querySpecIt = species.begin(); querySpecIt != species.end(); ++querySpecIt) {
	    /*
	      std::cout << "verifying inversion: " << *specIt 
	      << " " << *querySpecIt << std::endl;
	    */
	    VerifyInversion(query,
			    *specIt, *querySpecIt, seqName,
			    tempTableNames, maxTempTables, 
			    specStart[*specIt], specEnd[*specIt], 1, validation, seqMap);
	    
	    valInv->inversions.push_back(validation);
	    std::cout << " " << validation;
	  }
	  std::cout << std::endl;
	  (*binIt).second.push_back(valInv);
	}
      }
    }
    // Now the matrix for the inversion should be complete (every species is represented).
    //    PrintInversionListOrdered(species, (*binIt).second, std::cout);
    std::ofstream out;
    openck(outFileName, out);
    
    PrintBins(species, binnedInversions, out);
    out.close();
    ++binNumber;
  }
  std::ofstream out;
  openck(outFileName, out);

  PrintBins(species, binnedInversions, out);
  out.close();
}

