/***************************************************************************
 * Title:          joincontigs.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <vector>
#include <sstream>

#include "SeqReader.h"
#include "DNASequence.h"
#include "SeqUtils.h"

class CS {
public:
  ssize_t contig;
  ssize_t strand;
  CS& operator=(const CS &rhs) {
		if (this != &rhs) {
			contig = rhs.contig;
			strand = rhs.contig;
		}
    return *this;
  }
};

typedef std::vector<std::vector<CS>* > T_ChainsVect;

void ReadChains(std::ifstream &chainsFile, T_ChainsVect &chains);
void SkipWS(std::ifstream &in);
void BuildSuperContig(std::vector<DNASequence*> &contigs, std::vector<CS> chain, 
		      DNASequence &superContig);


void PrintUsage() {
  std::cout << "joincontigs contigs scaffolds outfile " << std::endl;
}

int main(int argc, char* argv[]) {

  std::vector<DNASequence*> contigs;

  std::ifstream in;
  if (argc < 2) {
    PrintUsage();
    exit(0);
  }
  in.open(argv[1]);

  SeqReader reader(&in);
  DNASequence *seq;
  while (reader.GetSeq(seq)) {
    contigs.push_back(seq);
  }

  in.close();
  in.clear();
  in.open(argv[2]);
  if (! in.good() ) {
    std::cout << " could not open " << argv[2] << std::endl;
    exit(0);
  }
  T_ChainsVect chains;
  ReadChains(in, chains);

  std::ofstream outfile;
  outfile.open(argv[3]);
  if (!outfile.good()) {
    std::cout << "could not open " << argv[3] << std::endl;
  }
  ssize_t i;
  DNASequence superContig;
  for (i = 0; i < chains.size(); i++) {
    BuildSuperContig(contigs, *chains[i], superContig);
    std::cout << "built contig of length: " << superContig.length << std::endl;
    //    superContig.PrintSeq(std::cout);
  }
  return 0;
}

void SkipWS(std::ifstream &in) {
  char c;
  while ((c = in.peek()) == ' ' || c == '\t') 
    in.get();
}

void ReadChains(std::ifstream &chainsFile, T_ChainsVect  &chains) {
  //UNUSED// char buf[2000];
  std::string stringBuf;
  std::vector< CS > *chain;
  CS chainStrand;
  while (!chainsFile.eof()) {
    chain = new std::vector<CS>;
    // Get the firs contig, always + dir
    chainsFile >> chainStrand.contig;
    chainStrand.strand = 1;
    chain->push_back(chainStrand);
    //UNUSED// ssize_t pos;
    char c;
    while ( (c = chainsFile.peek()) != EOF && c != '\n') {
      SkipWS(chainsFile);
      if (chainsFile.peek() == 'C') {
	chainStrand.strand = -1;
	chainsFile.get();
      }
      else {
	chainStrand.strand = 1;
      }
      chainsFile >> chainStrand.contig;
      chain->push_back(chainStrand);
      SkipWS(chainsFile);
    }
    chainsFile.get();
    std::cout << std::endl;
    chains.push_back(chain);
  }
}

void BuildSuperContig(std::vector<DNASequence*> &contigs, std::vector<CS> chain, 
		      DNASequence &superContig) {

  ssize_t i, j, pos;
  ssize_t contigLength = 0;
  ssize_t usedContigs = 0;
  DNASequence *contig;
  for (i = 0; i < chain.size(); i++) {
    if (chain[i].contig < contigs.size()) {
      contigLength += ((DNASequence*)contigs[chain[i].contig])->length;
      if (i < chain.size() - 1)
	contigLength += 100;
    }
  }
  superContig.Reset(contigLength);
  pos = 0;
  for (i = 0; i < chain.size(); i++) {
    if (chain[i].contig < contigs.size()) {
      contig = (DNASequence*) contigs[chain[i].contig];
      usedContigs++;
      if (chain[i].strand == 1) {
	for (j = 0; j < contig->length; j++) {
	  superContig.seq[pos] = contig->seq[j];
	  pos++;
	}
      }
      else {
	for (j = contig->length-1; j >= 0; j--) {
	  superContig.seq[pos] = NucRC(contig->seq[j]);
	  pos++;
	}
      }
      if (i < chain.size() - 1) {
	for (j = 0; j < 100; j++) {
	superContig.seq[pos] = 4;
	pos++;
	}
      }
    }
  }
  std::cout << "used: " << usedContigs << " contigs ";
}
