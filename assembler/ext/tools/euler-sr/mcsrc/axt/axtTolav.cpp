/***************************************************************************
 * Title:          axtTolav.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "axt/AxtReader.h"
#include "axt/AxtLib.h"
#include "lav/LAVFile.h"
#include "lav/LAVPrinter.h"

#include <iostream>
#include <vector>
#include <string>



void InitEnv(int argc, char* argv[], 
	     std::string &axtFileName, 
	     std::string &lavFileName,
	     ssize_t &refSeqLen,
	     ssize_t &qrySeqLen);

void PrintUsage();

void InitEnv(int argc, char* argv[], 
	     std::string &axtFileName,
	     std::string &lavFileName,
	     ssize_t &refSeqLen,
	     ssize_t &qrySeqLen) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "r:q:")) != EOF){
    switch(copt) {
    case 'r':
      refSeqLen = atoi(optarg);
      continue;
    case 'q':
      qrySeqLen = atoi(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify an axt file and lav file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  axtFileName = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify a lav file. " << std::endl;
    PrintUsage();
    exit(1);
  }
  lavFileName = argv[i];
}


void PrintUsage() {
  std::cout << "axtToLav.  Convert an axt file to a lav file. " << std::endl;
  std::cout << "usage: axtToLav [-r refseqlen] [-q qryseqlen] axtFileIn lavFileOutx " << std::endl;
}


int main(int argc, char *argv[]) {

  std::string axtFileName, lavFileName;
  ssize_t refSeqLen, qrySeqLen;
  refSeqLen = qrySeqLen = -1;
  InitEnv(argc, argv, axtFileName, lavFileName, refSeqLen, qrySeqLen);
  
  AxtEntries axtEntries;
  AxtReader::ReadAxtFile(axtFileName, axtEntries);
  //  PrintAxt(axtEntries, std::cout);
  LAVFile lavFile;
  AxtToLAV(axtEntries, lavFile, refSeqLen, qrySeqLen);

  std::ofstream out;
  out.open(lavFileName.c_str());
  if (!out.good()) {
    std::cout << "could not open " << lavFileName << std::endl;
    exit(0);
  }
  
  LAVPrinter::PrintLAVFile(lavFile, out);   
  out.close();
  return 0;
}
