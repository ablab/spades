/***************************************************************************
 * Title:          WindowCmp.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <ostream>
#include <stdio.h>

#include "DNASequence.h"
#include "SeqReader.h"
#include "alignutils.h"

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName);

void PrintUsage();

int main(int argc, char* argv[]) {

  std::string qrySeqName, refSeqName;
  DNASequence refSeq, qrySeq; 
  qrySeqName = "";
  refSeqName = "";
  InitEnv(argc, argv, refSeqName, qrySeqName);
  // Get input.
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);


  return 0;
}



void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &coordsFile) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "")) != EOF){
    switch(copt) {
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    PrintUsage();
    exit(1);
  }
  refSeqName = argv[i];
  i++;
  if (i >= argc) {
    PrintUsage();
    exit(1);
  }
  qrySeqName = argv[i];
}


void PrintUsage() {
  std::cout << "stub.  your description here. " << std::endl;
  std::cout << "usage: your usage here " << std::endl;
}
