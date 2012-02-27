/***************************************************************************
 * Title:          ValidateInversions.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <set>

// 3rd party
#include "mysql/mysql.h"

// Mine 
#include "net/NetDB.h"
#include "blocks/BlockDB.h"
#include "Inversion.h"
#include "InversionChars.h"
#include "InversionFile.h"
#include "InversionUtils.h"
#include "InversionDB.h"
#include "emboss/EmbossAlign.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "TupleLib.h"
#include "blocks/dbOrthoPos.h"
#include  "alignutils.h"

void CopySet(std::set<std::string> &setA, std::set<std::string> &setB) {
  std::set<std::string>::iterator setAIt;
  for (setAIt = setA.begin(); setAIt != setA.end(); ++setAIt) {
    setB.insert(*setAIt);
  }
}

void FindAllSpecies(InversionMap &invMap, StringSet &species);


ssize_t EnumerateInversions(SpeciesInversions &inversions, 
			std::vector<ssize_t> &invEnumerations,
			std::vector<ssize_t> &startPos,
			std::vector<ssize_t> &endPos,
			std::vector<std::string> &sequences) {

  // inversions is a list of inversions between this species and a
  // bunch of other species.  

  // Try to find out which inversions are overlapping to see how many
  // individual inversions this species particiaptes in.

  // invEnumerations is a vector of length \sum_i(n_i), where n_i is the number of species
  // that an inversion is found in.  

  SpeciesInversions::iterator invIt;

  ssize_t i, j, found;
  ssize_t invNumber = 0;
  ssize_t total = 0;
  ssize_t size = 0;
  invEnumerations.clear();
  startPos.clear();
  endPos.clear();

  for (invIt = inversions.begin(); invIt != inversions.end(); ++invIt) {
    size += (*invIt).second.size();
  }
  ssize_t notfound, totfound;
  notfound = 0;
  totfound = 0;
  for (invIt = inversions.begin(); invIt != inversions.end(); ++invIt) {
    found = 0;
    ssize_t start, end;
    for (i = 0; i < (*invIt).second.size();i++) {
      start = ((*invIt).second)[i]->tStart;
      end   = ((*invIt).second)[i]->tEnd;
      found = 0;
      ++total;
      for (j = 0; j < startPos.size() && ! found; j++) {
	// Look to see if this inversion overlaps any recorded inversions
	if ((start >= startPos[j] && start <= endPos[j]) or
	    (end >= startPos[j]   && end <= endPos[j]) or
	    (start <= startPos[j] && end >= endPos[j])) {
	  // Now compare the sizes of the two inversions to see if they
	  // are similar enough. 
	  if ((double(end - start) / double(endPos[j] - startPos[j])) > 0.4  and
	      (double(endPos[j] - startPos[j]) / double(end - start)) > 0.4) {
	    // This inversion is a valid overlap.
	    if (start < startPos[j])
	      startPos[j] = start;
	    if (end > endPos[j])
	      endPos[j] = end;
	  
	    // Annotate this inversion
	    invEnumerations.push_back(j);
	    found = 1;
	    totfound++;
	  }
	}
      }
      if (found == 0) {
	startPos.push_back(start);
	endPos.push_back(end);
	invEnumerations.push_back(invNumber);
	sequences.push_back((*invIt).second[0]->sequence);
	notfound++;
	invNumber++;
	found = 1;
      }
    }
  }
  return invNumber;
}

void GatherSpecies(SpeciesInversions &inversions, 
		   std::vector<ssize_t> &enumeration, 
		   ssize_t invEnum,
		   std::set<std::string> &with) {
  // Find out what species have an inversion, and what don't.  This is a 
  // rather inefficient way of finding them, but speed here doesn't matter.
  
  // enumeration has the inversion number for each species.  Look
  // through each inversion, and check to see if it's enumeration is 
  // invEnum.  If it is, add it to with, if not, without

  ssize_t i, j;
  ssize_t n = 0;
  SpeciesInversions::iterator invIt;
  std::vector<Inversion*> *invVectP;
  ssize_t size= 0;

  for (invIt = inversions.begin(); invIt != inversions.end(); invIt++) {
    size += (*invIt).second.size();
  }
  for (invIt = inversions.begin(); invIt != inversions.end(); invIt++) {
    invVectP = &((*invIt).second);
    for (i = 0; i < invVectP->size(); i++) {
      if (enumeration[n] == invEnum) {
	if (with.find((*invVectP)[i]->species) == with.end())
	  with.insert((*invVectP)[i]->species);
      }
      n++;
    }
  }
}

ssize_t verbose;

int main(int argc, char* argv[]) {

  std::string dbName;
  std::string tableName;
  std::string inversionFile;
  std::string validatedFileName;
  ssize_t maxNumTables = 1024;
  if (argc < 4) {
    std::cout << "usage: validateInversion database_name inversionFileName outputFileName [-nolocal] [-v] "
	      << std::endl;
    std::cout << "dbname and invname are as usual " << std::endl;
    std::cout << "validate inversions will take unknown inversions and decide if they " << std::endl;
    std::cout << "exist or not in all sequences.  Local alingment is performed in cases " << std::endl
	      << "where blastz alignments are not available, unless -nolocal is suplied" << std::endl;
    std::cout << "-v is for debugging " << std::endl;
    exit(1);
  }
  verbose = 0;
  dbName = argv[1];
  inversionFile = argv[2];
  validatedFileName = argv[3];
  
  int argi = 4;
  ssize_t doLocal = 1;
  ssize_t computeFullMat = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-v") == 0)
      ++verbose;
    else if (strcmp(argv[argi], "-nolocal") == 0)
      doLocal = 0;
    else if (strcmp(argv[argi], "-fullmat") == 0)
      computeFullMat = 1;
    ++argi;
  }
  std::ofstream validatedFile;
  openck(validatedFileName, validatedFile, std::ios::out);
  std::string scoreMatName;
  FloatMatrix scoreMat;

  char *home = getenv("HOME");
  scoreMatName = std::string(home) + "/projects/mcsrc/align/data/scoremat.txt";

  ReadScoreMatFile(scoreMatName, scoreMat);

  
  // Initialize everything.

  // Read the inversions
  InversionMap invMap;
  InversionMap::iterator invMapIt;
  StringSet species, with, without, unknown, gapped;
  std::vector<ssize_t> enumerations;
  std::vector<ssize_t> startPos, endPos;
  StringVector sequences;

  ReadInversionFile(inversionFile, invMap, species);

  // Connect to the database
  MYSQL *dbConn;
  ConnectToDB(dbConn, dbName);

  // Create the query and associate it with the databae
  DBQuery query(dbConn);


  std::string lavTableName, chainTableName, netTableName, tempTableName;

  StringSet tempTableNames;
  std::map<std::string, DNASequence *> seqMap;
  ssize_t enumerant;
  // build the list of species
  FindAllSpecies(invMap, species);

  StringSet::iterator sit;
  validatedFile << "species\tstart\tend\t";
  for (sit = species.begin(); sit != species.end(); ++sit) 
    validatedFile << "\t" << *sit;
  validatedFile << std::endl;
  // Invmap maps each reference specie to a list of species that it
  // has an inversion with, so invmap[human][chimp] is a list of
  // inversions between human and chimp.

  for (invMapIt = invMap.begin(); invMapIt != invMap.end(); ++invMapIt) {
    validatedFile << invMapIt->first << std::endl; 

    // For each reference species, find out how many inversions in the
    // query species overlap. Enumerate each non-overlapping inversion
    // by storing a value in enumerations, and record the start and
    // end position of each nonoverlapping inversion in the vectors
    // {start,end}Pos. 
    
    enumerations.clear();
    startPos.clear();
    endPos.clear();

    // Find out what inversions are at the same positions.
    ssize_t numInv;

    // enumerations contains the number of the inversion location,
    // where each inversion location is unique.
    numInv = EnumerateInversions(invMapIt->second, enumerations, startPos, endPos, sequences);
    ssize_t invStart, invEnd;

    // For each unique inversion, validate its presence on all query
    // sequences. 
    std::set<std::string>::iterator withoutIt, prevIt, removeIt;
    std::string refSeqName, qrySeqName;
    DNASequence* refSeq, *qrySeq;
      
    for (enumerant = 0; enumerant < numInv; enumerant++) {

      // Make it easier to reference the coordinates of the inversion
      invStart = startPos[enumerant];
      invEnd   = endPos[enumerant];

      // Reset the containers for inversion markers
      with.clear();
      without.clear();
      unknown.clear();
      gapped.clear();

      // Find what species share an inversion.  This finds all 0/1 (or +1/-1) relations
      // in the set of species.  
      GatherSpecies(invMapIt->second, enumerations, enumerant, with);

      if (verbose) {
	// Print with
	std::cout << invMapIt->first << " " << invStart << " " << invEnd << " with:";
	for (sit = with.begin(); sit != with.end(); ++sit) {
	  std::cout << " " << *sit;
	}
	std::cout << std::endl;
      }

      // Find what species don't share the inversion
      std::insert_iterator<std::set<std::string> > excluded(without, without.begin());  
      std::set_difference(species.begin(), species.end(),
			  with.begin(), with.end(),
			  excluded);
      if (verbose) {
	std::cout << "without: ";
	for (sit = without.begin(); sit != without.end(); ++sit) { 
	  std::cout << " " << *sit;
	}
	std::cout << std::endl;
      }

      // Now try to validate the inversions in the sequences that do
      // not have a recorded inversion.

      refSeqName = invMapIt->first + "." + sequences[enumerant] + ".fa";
      if (seqMap.find(refSeqName) == seqMap.end()) {
	refSeq = new DNASequence;
	SeqReader::GetSeq(refSeqName, *refSeq, SeqReader::noConvert);
	seqMap[refSeqName] = refSeq;
      }
      else {
	refSeq = seqMap[refSeqName];
      }

      std::string spec;
      ssize_t validation;
      for (withoutIt = without.begin(); withoutIt != without.end(); ) {
	// Configure names of tables and sequence files.
	VerifyInversion(query,
			invMapIt->first, *withoutIt, sequences[enumerant],
			tempTableNames, maxNumTables, 
			startPos[enumerant], endPos[enumerant], doLocal, 
			validation, seqMap,
			scoreMat);

	if (validation == WITH_INVERSION) {
	  spec = *withoutIt;
	  prevIt = withoutIt;
	  ++withoutIt;
	  without.erase(prevIt);
	  with.insert(spec);
	}
	else if (validation == WITHOUT_INVERSION) {
	  ++withoutIt;
	}
	else if (validation == UNKNOWN) {
	  spec = *withoutIt;
	  prevIt = withoutIt;
	  ++withoutIt;
	  without.erase(prevIt);
	  unknown.insert(spec);
	}
	else if (validation == DELETED) {
	  spec = *withoutIt;
	  prevIt = withoutIt;
	  ++withoutIt;
	  without.erase(prevIt);
	  gapped.insert(spec);
	}
	else {
	  std::cout << "Did not assign inversion a category " << std::endl;
	}
      }

      // Now print a summary of the alignments for this species
      validatedFile << invMapIt->first << "\t" << invStart << "\t" << invEnd;
      for (sit = species.begin(); sit != species.end(); ++sit) {
	validatedFile << "\t";
	if (with.find(*sit) != with.end()) 
	  validatedFile << WITH_INVERSION;
	else if (without.find(*sit) != without.end()) 
	  validatedFile << WITHOUT_INVERSION;
	else if (unknown.find(*sit) != unknown.end())
	  validatedFile << UNKNOWN;
	else if (gapped.find(*sit) != gapped.end())
	  validatedFile << DELETED;
	else {
	  validatedFile << "NC" << std::endl;
	}
      }
      validatedFile << std::endl;
    }
  }
  return 0;
}

void FindAllSpecies(InversionMap &invMap, StringSet &species) {
  InversionMap::iterator invMapIt;
  SpeciesInversions::iterator specIt;

  for (invMapIt = invMap.begin(); invMapIt != invMap.end(); ++invMapIt) {
    if (species.find(invMapIt->first) == species.end()) {
      species.insert(invMapIt->first);
    }
    for (specIt = invMapIt->second.begin(); specIt != invMapIt->second.end(); ++specIt) {
      if (species.find(specIt->first) == species.end()) {
	species.insert(specIt->first);
      }
    }
  }
}

