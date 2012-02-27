/***************************************************************************
 * Title:          stripsToFasta.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <string>
#include <stdlib.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sstream>

#include "StripGen.h"
#include "SeqReader.h"
#include "DNASequence.h"
#include "TupleLib.h"
#include "SeqUtils.h"

void PrintUsage();

void initenv(int argc, char *argv[], 
	     std::string &refFile,
	     std::string &qryFile,
	     std::string &stripFile,
	     std::string &refOutFile,
	     std::string &qryOutFile);


int main(int argc, char* argv[]) {
  std::string refFileName, qryFileName, stripFileName, refOutName, qryOutName;
  std::string indexA("@");
  std::string indexB("!");
  initenv(argc, argv, refFileName, qryFileName, stripFileName, refOutName, qryOutName);

  T_Strips strips;
  
  std::ifstream stripFile;
  
  stripFile.open(stripFileName.c_str());
  if ( ! stripFile.good() ) {
    std::cout << "Could not open infile " << std::endl;
    exit(0);
  }
  ssize_t numRead;
  GetStrips(stripFile, strips, numRead, 0);

  std::ifstream refFile, qryFile;

  refFile.open(refFileName.c_str());
  qryFile.open(qryFileName.c_str());
  
  DNASequence refSeq, qrySeq;

  SeqReader::GetSeq(refFile, refSeq,1);
  SeqReader::GetSeq(qryFile, qrySeq,1);
  
  refFile.close();
  qryFile.close();
  T_Strips::iterator it,end;
  end = strips.end();
  ssize_t step, qryPos, refPos;
  ssize_t startQry, endQry, startRef, endRef;
  DNASequence qryStrip, refStrip;
  DNASequence qryRCStrip;
  std::ostringstream seqNameStream;
  std::string seqName;
  ssize_t stripNumber, stripLength;
  char refDir, qryDir;

  std::ofstream refOut, qryOut;
  refOut.open(refOutName.c_str());
  qryOut.open(qryOutName.c_str());

  if (! refOut.good() ) {
    std::cout << "could not open " << refOutName << std::endl;
    exit(0);
  }
  if (! qryOut.good() ) {
    std::cout << "could not open " << qryOutName << std::endl;
    exit(0);
  }

  for (it = strips.begin(); it != strips.end(); ++it) {
    // Print the reference strip.  It is always assumed this is in 
    // the (+) dir.
    refStrip.seq = &refSeq.seq[it->startRefPos];
    refStrip.length = it->endRefPos - it->startRefPos + 1;
    seqName = "";
    seqNameStream.str(seqName);
    if (it->endRefPos - it->startRefPos > 0) 
      refDir = '+';
    else
      refDir = '-';
    if (it->endQryPos - it->startQryPos > 0)
      qryDir = '+';
    else
      qryDir = '-';
    
    seqNameStream << it->startQryEnum << " " << refStrip.length << " " << refDir << " " << qryDir;
    refStrip.StoreName((char*) seqNameStream.str().c_str());
    refStrip.PrintSeq(refOut);
    
    // Print the query sequence.  It's possible that this sequence
    // is a reverse complement, so it's a bit more of a pain to print.
    if (it->endQryPos - it->startQryPos > 0) {
      // Print in (+) dir
      qryStrip.seq = &qrySeq.seq[it->startQryPos];
      qryStrip.length = it->endQryPos - it->startQryPos + 1;
      seqName = "";
      seqNameStream.str(seqName);
      seqNameStream << it->startQryEnum << " " << qryStrip.length << " " << qryDir << " " << refDir;
      qryStrip.StoreName((char*) seqNameStream.str().c_str());
      qryStrip.PrintSeq(qryOut);
    }
    else {
      // Print in (-) dir
      qryStrip.seq = &qrySeq.seq[it->startQryPos];
      qryStrip.length = it->startQryPos - it->endQryPos + 1;
      MakeRC(qryStrip, qryRCStrip);
      seqName = "";
      seqNameStream.str(seqName);
      seqNameStream << it->startQryEnum << " " << qryRCStrip.length << " " << qryDir << " " << refDir;
      qryRCStrip.StoreName((char*)seqNameStream.str().c_str());
      qryRCStrip.PrintSeq(qryOut);
    }
  }

  qryOut.close();
  refOut.close();
  stripFile.close();
}

void initenv(int argc, char *argv[], 
	     std::string &refFile,
	     std::string &qryFile,
	     std::string &stripFile,
	     std::string &refOutFile,
	     std::string &qryOutFile) {
  ssize_t copt;
  if (argc != 5) {
    PrintUsage();
    exit(1);
  }
  refFile = argv[1];
  qryFile = argv[2];
  stripFile = argv[3];
  refOutFile = argv[4];
  qryOutFile = argv[5];
}
void PrintUsage() {
  std::cout << "stripsToFasta, a program to generate fasta entries for the strips in an alignment" << std::endl;
  std::cout << "usage:  stripsToFasta ref.fasta qry.fasta enumFile refOut.fasta qryOut.fasta" << std::endl; 
}







