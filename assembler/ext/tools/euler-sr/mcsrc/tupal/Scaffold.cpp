/***************************************************************************
 * Title:          Scaffold.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "TupalLib.h"
#include <iostream>

void PrintUsage();

void InitEnv(int argc, char* argv[], 
	     std::string &refSeqName, 
	     std::string &qrySeqName,
	     std::string &scaffoldName, 
	     T_Params &params, T_ParallelParams parInfo) {
  ssize_t copt;
  ssize_t i;
  std::string inpfile;
  while ( (copt=getopt(argc, argv, "h:w:r:")) != EOF){
    switch(copt) {
    case 'h':
      params.hashLength = atoi(optarg);
      continue;
    case 'w':
      params.wordLength = atoi(optarg);
      continue;
    case 'r':
      ReadParamFile(optarg, params, parInfo);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a reference seq and query seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  refSeqName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a query seq. " << std::endl;
    PrintUsage();
    exit(1);
  }
  qrySeqName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a scaffold. " << std::endl;
    PrintUsage();
    exit(1);
  }
  scaffoldName = argv[i];
}

void PrintUsage() {
  std::cout << "scaffold, a program to generate a subset of conserved, unique matches " << std::endl
	    << "between two sequences such that the set is a longest in(de)creasing subset " << std::endl
	    << "from all conserved unique matches. " << std::endl << std::endl;
  std::cout << "usage: [-h hashLength ] [-w wordLength] [-r paramfile] refSeq qrySeq scaffoldName  " << std::endl;
}


int main(int argc, char *argv[]) {
  std::string refSeqName, qrySeqName, scaffoldName;
  T_Params params;
  T_ParallelParams parInfo;
  DNASequence refSeq, qrySeq;
  T_Alignment *alignment;
  T_Strips* scaffold;
  InitEnv(argc, argv, refSeqName, qrySeqName, scaffoldName, params, parInfo);
  
  SeqReader::GetSeq(refSeqName, refSeq);
  SeqReader::GetSeq(qrySeqName, qrySeq);


  BuildScaffold(refSeq, qrySeq, alignment, scaffold, params, parInfo,"");

  std::ofstream out;
  out.open(scaffoldName.c_str());
  if (!out.good() ) {
    std::cout << "Could not open scaffoldling output. " << std::endl;
    exit(1);
  }
  PrintEnumerations(*scaffold, out);
  
  if (alignment != NULL)
    FreeAlignments(alignment);
  alignment = NULL;
  delete[] refSeq.seq;
  delete[] qrySeq.seq;
  
  
  return 0;
}
  
