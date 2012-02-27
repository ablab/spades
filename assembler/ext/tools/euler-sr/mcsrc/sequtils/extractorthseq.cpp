/***************************************************************************
 * Title:          extractorthseq.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <fstream>

#include <unistd.h>
#include <assert.h>
#include <sstream>
#include "utils.h"
#include "SeqReader.h"
#include "DNASequence.h"



void initenv(int argc, char *argv[], 
	     std::string &seqADir, std::string &seqBDir,
	     std::string &positionfile,
	     std::string &outfile);

void printusage();


class Cpos {
public:
  std::string fileA, fileB;
  ssize_t startA, endA, startB, endB;
  Cpos & operator=(const Cpos &rhs) {
		if (this != &rhs) {
			fileA  = rhs.fileA;
			fileB  = rhs.fileB;
			startA = rhs.startA;
			endA   = rhs.endA;
			startB = rhs.startB;
			endB   = rhs.endB;
		}
    return *this;
  }
};

typedef std::vector<Cpos > PosVector;

int main(int argc, char *argv[]) {
  std::ofstream outFile;
  std::string seqFileName, positionFileName, outBase;
  std::string seqADir, seqBDir;

  initenv(argc, argv, seqADir, seqBDir, positionFileName, outBase);
  PosVector positions;
  std::ifstream positionFile;
  openck(positionFileName, positionFile);

  while (positionFile) {
    Cpos pos;
    positionFile >> pos.fileA >> pos.startA >> pos.endA 
		 >> pos.fileB >> pos.startB >> pos.endB;
    positions.push_back(pos);
  }

  ssize_t p;
  std::string seqAName, seqBName, seqAFileName, seqBFileName;
  DNASequence seqA, seqB, subseqA, subseqB;
  if (positions.size() == 0) 
    return 0;

  seqAName = "";
  seqBName = "";
  subseqA._ascii = 1;
  subseqB._ascii = 1;
  for (p = 0; p < positions.size(); p++ ) {
    std::cout << "extracting: " <<  positions[p].fileA << " " <<  positions[p].startA << " " << positions[p].endA  << " " <<  positions[p].fileB << " " <<  positions[p].startB << " " <<  positions[p].endB << std::endl;
    if (positions[p].fileA != seqAName) {
      seqAName = positions[p].fileA;
      std::cout << "reading " << seqAName << std::endl;
      seqAFileName = seqADir + "/" + seqAName + ".fa";
      SeqReader::GetSeq(seqAFileName, seqA, SeqReader::noConvert);
      std::cout << " done " << std::endl;
    }
    if (positions[p].fileB != seqBName) {
      seqBName = positions[p].fileB;
      std::cout << "reading " << seqBName << std::endl;
      seqBFileName = seqBDir + "/" + seqBName + ".fa";
      SeqReader::GetSeq(seqBFileName, seqB, SeqReader::noConvert);
      std::cout << " done " << std::endl;
    }
    
    subseqA.seq = &(seqA.seq[positions[p].startA]);
    subseqA.length = positions[p].endA - positions[p].startA;
    std::cout << "alength: " << subseqA.length << std::endl;
    subseqB.seq = &(seqB.seq[positions[p].startB]);
    subseqB.length = positions[p].endB - positions[p].startB;
    std::cout << "blength: " << subseqB.length << std::endl;

    std::ofstream outA, outB;
    std::stringstream outAStrm, outBStrm;
    std::string outAName, outBName;
    outAStrm << outBase << "." << positions[p].fileA << "." << positions[p].startA << "." << positions[p].endA << ".fasta";
    outBStrm << outBase << "." << positions[p].fileB << "." << positions[p].startB << "." << positions[p].endB << ".fasta";

    outAName = outAStrm.str();
    outBName = outBStrm.str();

    subseqA.namestr = outAName;
    openck(outAName, outA, std::ios::out);
    subseqA.PrintSeq(outA);
    outA.close();
    outA.clear();

    subseqB.namestr = outBName;
    openck(outBName, outB, std::ios::out);
    subseqB.PrintSeq(outB);
    outB.close();
    outB.clear();
  }
  return 0;
}


void initenv(int argc, char *argv[], 
	     std::string &seqADir, std::string &seqBDir,
	     std::string &positionFile,
	     std::string &outbase){
  //UNUSED// ssize_t copt;
  if (argc != 5) { 
    printusage();
    exit(0);
  }
  seqADir = argv[1];
  seqBDir = argv[2];
  positionFile = argv[3];
  outbase = argv[4];
}

void printusage() {
  std::cout << "usage:  extractorthseq seqadir seqbdir posfile outbase " << std::endl;
}
