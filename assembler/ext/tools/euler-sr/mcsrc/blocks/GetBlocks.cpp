/***************************************************************************
 * Title:          GetBlocks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <DNASequence.h>
#include <StripGen.h>


void initenv(int argc, char *argv[], 
	     std::string &enumFile,
	     std::string &outputFile,
	     ssize_t &gapLength);
void PrintUsage();

int main(int argc, char* argv[]) {

  std::string stripsFileName;
  std::string blocksFileName;
  ssize_t gapLength = 10;
  T_Strips strips;
  initenv(argc, argv, stripsFileName, blocksFileName, gapLength);

  std::ifstream stripsFile;
  stripsFile.open(stripsFileName.c_str());
  if (!stripsFile.good()) {
    std::cout << "Could not open " << stripsFileName << std::endl;
    exit(0);
  }
  ssize_t numStrips;
  GetStrips(stripsFile, strips, numStrips, 0);

  std::cout << "got " << numStrips << " strips " << std::endl;
}


void initenv(int argc, char *argv[], 
	     std::string &stripsFile,
	     std::string &outputFile,
	     ssize_t &gapLength) 
{
  ssize_t copt;
  ssize_t i;
  std::string inpfile;
  while ( (copt=getopt(argc, argv, "g:")) != EOF){
    switch(copt) {
    case 'g':
      gapLength = atoi(optarg);
      continue;
    }
  }
  i = optind;
  if (i >= argc) {
    std::cout << "You must specify a strips file. " << std::endl;
    PrintUsage();
    exit(0);
  }
  stripsFile = argv[i];
  i++;
  if (i >= argc) {
    std::cout << "You must specify an output blocks file. " << std::endl;
    PrintUsage();
    exit(0);
  }
}


void PrintUsage() {
  std::cout << "getblocks.  A program to find the contiguous blocks from " << std::endl
	    << " two sequences given an enumeration of the aligned locations. " << std::endl
	    << "A contiguous sequence is a sequence of positions mapped " << std::endl
	    << "with up to gapLength gap allowed" << std::endl;
  std::cout << "usage: getblocks [-g gapLength] strips_in blocks_out " << std::endl;
}


    
       
