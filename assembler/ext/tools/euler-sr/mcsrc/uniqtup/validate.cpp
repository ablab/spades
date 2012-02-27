/***************************************************************************
 * Title:          validate.cpp 
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
#include <fstream>
#include <stdlib.h>

#include "SeqReader.h"
#include "TupleLib.h"
#include "DNASequence.h"
#include "utils.h"


void printusage();
void initenv(int argc, char *argv[], 
	     std::vector< std::string> &probes,
	     std::vector< std::string> &probeFiles,
	     std::vector<std::string> &infile,
	     std::string &outfile);

void ReadProbeFile(std::string probeFile, std::vector<std::string> &probes) {

  std::ifstream in;
  openck(probeFile, in);
  std::string probe;
  ssize_t count; 
  ssize_t others;
  while (in) {
    if ( ! ( in >> count >> probe >> count) )
      break;
    else
      probes.push_back(probe);
  }
  in.close();
}

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
		  std::ostream &out) { 
  ssize_t i;
  for (i = 0; i < length; i++) {
    if (mult[i] != 0) {
      out  << mult[i] << " " ;
      std::cout << "printing tuple: " << tuples[i] << std::endl;
      printTuple(tuples[i], tupleLength, out);
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
  std::vector<std::string> probes, probesRC, probeFileNames;

  initenv(argc, argv, 
	  probes, probeFileNames, infileNames, outfileExt);
  
  for (i = 0; i < probeFileNames.size(); i++) {
    ReadProbeFile(probeFileNames[i], probes);
  }
  DNASequence rc;
  ssize_t **probeCounts;
  for (i = 0; i < probes.size(); i++) {
    MakeRC((char*) probes[i].c_str(), probes[i].size(), rc.seq);
    rc.length = probes[i].size();
    std::cout << "adding probe of len: " << rc.length << std::endl;
    probes.push_back(rc.seq);
  }
  probeCounts = new ssize_t*[probes.size()];

  for (i = 0; i < probes.size(); i++) {
    probeCounts[i] = new ssize_t[infileNames.size()];
    for (j = 0;j < infileNames.size(); j++)
      probeCounts[i][j] = 0;
  }

  
  DNASequence seq;
  ssize_t p, pos;
  ssize_t probeLen;
  probeLen = probes[0].size();
  for (i  = 0; i < infileNames.size(); i++) {
    std::cout << infileNames[i] << std::endl;
    SeqReader::GetSeq(infileNames[i], seq, SeqReader::noConvert);
    for (pos = 0; pos < seq.length; pos++) {
      for (p = 0; p < probes.size(); p++) {
	if (strncmp(probes[p].c_str(), &seq.seq[pos], probeLen) == 0)
	  probeCounts[p][i]++;
      }
    }
  }

  for (i = 0; i < infileNames.size(); i++)
    std::cout << infileNames[i] << ",";
  std::cout << std::endl;
  for (i = 0; i < probes.size() / 2; i++) {
    std::cout << probes[i] << ",";
    for (j = 0; j < infileNames.size(); j++)
      std::cout << probes[i][j] + probes[i+probes.size()/2] << ",";
    std::cout << std::endl;
  }
  return 0;
}

void initenv(int argc, char *argv[], 
	     std::vector<std::string> &probes,
	     std::vector<std::string> &probeFiles,
	     std::vector<std::string> &infileNames,
	     std::string &outfile) {
  ssize_t i;
  ssize_t copt;
  std::string infile;
  // do my own parsing of the command line
  i = 1;
  ssize_t searchProbes = 1;
  while (i < argc) {
    if (strcmp(argv[i], "-f") == 0 || strcmp(argv[i], "-p") == 0) 
      break;
    else
      probes.push_back(argv[i]);
    i++;
  }
  if (strcmp(argv[i], "-p") == 0) {
    i++;
    while (i < argc) {
      if (strcmp(argv[i], "-f") == 0)
	break;
      else
	probeFiles.push_back(argv[i]);
      i++;
    }
  }
  if (i >= argc) {
    std::cout << "you must enter -f filenames " << std::endl;
    exit(0);
  }
  while (i < argc ) {
    if (strcmp(argv[i], "-o") == 0)
      break;
    infileNames.push_back(argv[i]);
    i++;
  }
  if (i < argc) {
    outfile = argv[i];
  }
}

void printusage() {
  std::cout << "usage:  findunique  tuple1 ... -p probefile1 -i tuple_count1 [-i tuplecount2 ...] -t tuplelength " << std::endl;  
  std::cout << "-i probes_file1 ... : files to look for probes in " << std::endl;
  std::cout << "-i input sequence      : sequences to find unique probes in." << std::endl; 
  std::cout << "-o output              : output extension." << std::endl; 
  std::cout << "-t tuplelength         : length of the probes to search for (20)." << std::endl; 
}
