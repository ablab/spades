/***************************************************************************
 * Title:          findunique.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <ios>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <ostream>
#include <stdlib.h>

#include "SeqReader.h"
#include "TupleLib.h"
#include "DNASequence.h"


void printusage();
void initenv(int argc, char *argv[], 
	     ssize_t& tupleLen, 
	     ssize_t& maxMultiplicity,
	     ssize_t& minQryMultiplicity,
	     std::vector<std::string> &infile,
	     std::string &outfile);


void GetMin(Tuple **tuples, ssize_t numarrays, ssize_t *pos, ssize_t &min) {
  // Find the column of tuples that has the minimum value
  // where pos is the position in each column to look.
  ssize_t i;
  min = tuples[0][pos[0]];
  ssize_t minpos = 0;
  for (i = 1; i < numarrays; i++) {
    if (tuples[i][pos[i]] < min) {
      min = tuples[i][pos[i]];
      minpos = i;
    }
  }
}

void GetEqualIndices(Tuple **tuples, ssize_t numarrays, 
		     ssize_t *pos, Tuple value, ssize_t *indices) {
  // Find all columns of tuples where the value at pos
  // is equal to the parameter 'value'.
  ssize_t i;
  for (i = 0; i < numarrays; i++) {
    if (tuples[i][pos[i]] == value)
      indices[i] = 1;
    else
      indices[i] = 0;
  }
}

void StoreUnique(Tuple **uniqueTuples, ssize_t **multiplicities,
		 ssize_t* pos, // all positions of unique tuples
		 ssize_t index, // the set of tuples to add to
		 Tuple value, ssize_t mult // the tuple and mult to add
		 ) {
  uniqueTuples[index][pos[index]] = value;
  multiplicities[index][pos[index]] = mult;
  pos[index]++;
}

void ReadSeq(std::ifstream &in, Tuple *&tuples, ssize_t *&mult, ssize_t& length) {
  ssize_t i;
  in.read((char*)&length, sizeof(ssize_t));
  tuples = new Tuple[length];
  mult   = new ssize_t[length];
  for (i = 0; i < length; i++ ) {
    in.read((char*)&mult[i], sizeof(ssize_t));
    in.read((char*)&tuples[i], sizeof(Tuple));
  }
}


void MaskUnique(Tuple* qryTuples, ssize_t *qryMult, ssize_t qryLength,
		Tuple* refTuples, ssize_t *refMult, ssize_t refLength,
		ssize_t minReferenceMult, ssize_t minMultiplicity) {
  // mask off all non-unique tuples in the query list.
  // masking is done by setting the multiplicity of a tuple to 0.
  ssize_t qryPos, refPos;
  qryPos = refPos = 0;
  ssize_t lastMasked;
  ssize_t numUnmasked = 0;
  while (refPos < refLength && qryPos < qryLength) {
    // advance past any masked query positions.

    while(qryPos < qryLength && 
	  (qryMult[qryPos] == 0 || qryMult[qryPos] < minMultiplicity)) {
      qryMult[qryPos] = 0; // just in case this position wasn't masked before (it had too low mult
      qryPos++;
      lastMasked = qryPos;
    }
    if (qryTuples[qryPos] == refTuples[refPos]) {
      if(qryMult[qryPos] < minMultiplicity) {
	std::cout << "query mult is : " <<  qryMult[qryPos] 
		  << " but min was: " << minMultiplicity << std::endl;
	std::cout << "the last masked position was: " << lastMasked 
		  << std::endl;
      }
      if (refMult[refPos] >= minReferenceMult) {
	// It is ok to keep a tuple if the multiplicity of the tuple
	// in the reference genome is below a certain point.  The
	// converse is that when there is above a certain threshold in
	// the reference genome, the tuple should be masked.
	qryMult[qryPos] = 0;
      }
      else {
	numUnmasked++;
      }
      qryPos++;
      refPos++;
    } else if (qryTuples[qryPos] > refTuples[refPos]) {
      refPos++;
    }
    else if (qryTuples[qryPos] < refTuples[refPos]) {
      numUnmasked++;
      qryPos++;
    }
  }
  std::cout << "  ended with: " << numUnmasked << " positions " <<std::endl;
}

void PrintNonZero(Tuple *tuples, ssize_t* mult, ssize_t length, ssize_t tupleLength,
		  std::ofstream &out) { 
  ssize_t i;
  for (i = 0; i < length; i++) {
    if (mult[i] != 0) {
      out  << mult[i] << " " ;
      std::cout << "printing tuple: " << tuples[i] << std::endl;
      Tuple t;
      t = 12345;
      printTuple(13, 12, out);
      //      printTuple(tuples[i], tupleLength, std::cout);
      out << std::endl ;
    }
  }
}
      
int main(int argc, char* argv[]) {


  std::vector<std::string> infileNames;
  std::string  outfileName;
  ssize_t tupleLen;
  ssize_t minRefMultiplicity;
  ssize_t i, j;
  ssize_t numSeq;
  ssize_t totalTuples;
  std::string outfileExt;
  ssize_t minQryMultiplicity;

  outfileExt = ".unique";
  minRefMultiplicity = 0;
  minQryMultiplicity = 0;
  tupleLen = -1;
    
  initenv(argc, argv, tupleLen, 
	  minRefMultiplicity, minQryMultiplicity,
	  infileNames, outfileExt);
  if (tupleLen == -1) {
    std::cout << "must specify a tuple length " << std::endl;
    exit(1);
  }

  std::ifstream qryIn, refIn;
  std::ofstream uniqueOut;
  std::string uniqueOutName;
  if (infileNames.size() == 0) {
    printusage();
    exit(0);
  }
  numSeq = infileNames.size();

  Tuple *qryTuples, *refTuples;
  ssize_t  *qryMult, *refMult;
  ssize_t    qryLen,   refLen;
    
  for (i = 0; i < infileNames.size(); i++) {
    qryIn.open(infileNames[i].c_str(), std::ios::binary);
    if (!qryIn.good()) {
      std::cout << " could not open " << infileNames[i] << std::endl;
      exit(1);
    }
    ReadSeq(qryIn, qryTuples, qryMult, qryLen);
    std::cout << infileNames[i] << ": " << std::endl;
    for (j = 0; j < infileNames.size(); j++) {
      if (i != j) {
	refIn.open(infileNames[j].c_str(), std::ios::binary);
	ReadSeq(refIn, refTuples, refMult, refLen);
	//	std::cout << " vs " <<  infileNames[j];
	MaskUnique(qryTuples, qryMult, qryLen,
		   refTuples, refMult, refLen,
		   minRefMultiplicity, minQryMultiplicity);

	delete refTuples;
	delete refMult;
	refIn.close();
	refIn.clear();
      }
    }
    qryIn.close();
    qryIn.clear();
    uniqueOutName = infileNames[i] + "." + outfileExt;
    uniqueOut.open(uniqueOutName.c_str());
    if (! uniqueOut.good() ) {
      std::cout << "Could not open " << uniqueOutName << std::endl;
      exit(0);
    }
    PrintNonZero(qryTuples, qryMult, qryLen, tupleLen, uniqueOut );
    uniqueOut.close();
    uniqueOut.clear();
    delete qryTuples;
    delete qryMult;
  }

  return 0;
}

    
ssize_t GetNext(Tuple* list, Tuple cur, ssize_t &pos) {
  ssize_t advanced = 0;
  while (list[pos] == cur) {
    advanced++;
    pos++;
  }  
  return list[pos];
}   

void initenv(int argc, char *argv[], 
	     ssize_t& tupleLen, 
	     ssize_t& minRefMultiplicity,
	     ssize_t& minQryMultiplicity,
	     std::vector<std::string> &infileNames,
	     std::string &outfile) {
  ssize_t copt;
  std::string infile;
  while ( (copt=getopt(argc, argv, "i:o:t:M:m:")) != EOF) {
    switch(copt) {
    case 'o':
      outfile = optarg;
      continue;
    case 't':
      sscanf(optarg,"%d", &tupleLen);
      continue;
    case 'M':
      sscanf(optarg, "%d", &minRefMultiplicity);
      continue;
    case 'm':
      sscanf(optarg, "%d", &minQryMultiplicity);
      continue;
    default:
      printusage();
      exit(1);
    }
  }
  ssize_t i;
  i = optind;
  while (i < argc) {
    infileNames.push_back(argv[i]);
    i++;
  }
  if (infileNames.size() == 0 or outfile == "") {
    printusage();
    exit(1);
  }
}

void printusage() {
  std::cout << "usage:  findunique  -i tuple_count1 [-i tuplecount2 ...] -t tuplelength " << std::endl;  
  std::cout << "-i input sequence      : sequences to find unique probes in." << std::endl; 
  std::cout << "-o output              : output extension." << std::endl; 
  std::cout << "-t tuplelength         : length of the probes to search for (20)." << std::endl; 
}
