/***************************************************************************
 * Title:          simulator.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>


void InitEnv(int argc, char* argv[], 
	     ssize_t &seqLen,
	     double &pctGC,
	     std::string &outFileName,
	     ssize_t lineLength=60);


int main(int argc, char* argv[]) {



  return 0;
}

void InitEnv(int argc, char* argv[], 
	     ssize_t &seqLen,
	     double &pctGC,
	     std::string &outFileName,
	     ssize_t lineLength=60) {

  while ( (copt=getopt(argc, argv, "l:g:")) != EOF) {
    switch(copt) {
    case 'l':
      lineLength = atoi(optarg);
      continue;
    case 'g':
      pctGC = atof(optarg);
      continue;
    default:
      PrintUsage();
      exit(1);
    }
  }
}

void PrintUsage() {
  std::cout << "simulator: create random sequences with gc bias " << std::endl;
  std::cout << "usage: "
}


