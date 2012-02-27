/***************************************************************************
 * Title:          sorttuples.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  02/28/2010
 *
 * Copyright (c) 2007-2010 The Regents of the University of California
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

// TODO: may need to update file format to have tupleSize at start, as .spect does

void printusage();
void initenv(int argc, char *argv[], 
	     ssize_t& probelen, 
	     ssize_t& minCount,
	     std::string &infile,
	     std::string &outfile);

ssize_t compare(const void* ap, const void* bp) {
  Tuple a, b;
  a = (Tuple) *((Tuple*)ap);
  b = (Tuple) *((Tuple*)bp);
  if (a == b) return 0;
  if (a < b) return -1;
  return 1;
}

int main(int argc, char* argv[]) {

  
  std::string infileName, outfileName;
  ssize_t tupleLen;
  ssize_t minCount;
  infileName  = "";
  outfileName = "";
  tupleLen = 20;
  minCount = 5;
  initenv(argc, argv, tupleLen,minCount, infileName, outfileName);
  std::ifstream in;

  in.open(infileName.c_str());
  if (!in.good()) {
    std::cout << "could not open " << infileName << std::endl;
    exit(1);
  }

  SeqReader reader(&in);
  std::vector<DNASequence*> sequences;
  ssize_t totalTuples = 0;
  DNASequence *seq;
  DNASequence *rev;
  char* revSeq;
  ssize_t len;

  while (reader.GetSeq(seq)) { 
    sequences.push_back(seq);
    totalTuples += (seq->length - tupleLen + 1);
    rev = new DNASequence;
    MakeRC(*seq, *rev);
    totalTuples += (rev->length - tupleLen + 1);
    sequences.push_back(rev);
  }

  std::cout << "done reading the seq " << totalTuples << std::endl;
  Tuple* tuples = new Tuple[totalTuples];
  //  Tuple* sortedTuples = new Tuple[totalTuples];
  ssize_t tupleIndex = 0;
  ssize_t i,j;
  Tuple tuple;
  for (i = 0; i < sequences.size(); i++) {
    for (j = 0; j < sequences[i]->length - tupleLen + 1; j++) {
      tuple = trans_seq(&(sequences[i]->seq[j]), tupleLen);
      tuples[tupleIndex] = tuple;
      tupleIndex++;
    }
  }
  qsort(tuples, totalTuples, sizeof(Tuple), compare);
  
  std::ofstream output;
  output.open(outfileName.c_str(), std::ios::binary);
  ssize_t unique = 0;
  output.write((char*) &unique, sizeof(ssize_t));
  //  output.write((char*)tuples, sizeof(Tuple)*totalTuples);
  ssize_t pos = 0;
  ssize_t advanced;
  Tuple cur;
  char mult;
  while (pos < totalTuples) {
    advanced = 0;
    cur = tuples[pos];
    mult = 1;
    while (pos < totalTuples && tuples[pos] == cur) {
      advanced++;
      pos++;
      mult++;
    }
    if (pos <= totalTuples) {
      /*
	printTuple(tuples[pos-1], tupleLen, std::cout);
	std::cout << std::endl;
      */
      if (mult > minCount) {
	unique++;
	output.write((char*) &advanced, sizeof(ssize_t));
	output.write((char*) &tuples[pos-1], sizeof(Tuple));
      }
    }
  }
  output.seekp(0);
  output.write((char*) &unique, sizeof(ssize_t));
  std::cout << unique << " unique out of " << totalTuples << std::endl;

  return 0;
}

void initenv(int argc, char *argv[], 
	     ssize_t& probelen, 
	     ssize_t& minCount,
	     std::string &infile,
	     std::string &outfile) {
  ssize_t copt;
  while ( (copt=getopt(argc, argv, "i:o:t:m:")) != EOF) {
    switch(copt) {
    case 'm':
      minCount = atoi(optarg);
      continue;
    case 'i':
      infile = optarg;
      continue;
    case 'o':
      outfile = optarg;
      continue;
    case 't':
      probelen = atoi(optarg);
      continue;
    default:
      printusage();
      exit(1);
    }
  }
  if (infile == "" or outfile == "") {
    printusage();
    exit(1);
  }
}

void printusage() {
  std::cout << "usage:  sorttuples  -i seq1 -t tuplelength " << std::endl;  
  std::cout << "-i input sequence      : sequences to find unique probes in." << std::endl; 
  std::cout << "-o output              : output file." << std::endl; 
  std::cout << "-t tuplelength         : length of the tuple to search for (20)." << std::endl; 
}
