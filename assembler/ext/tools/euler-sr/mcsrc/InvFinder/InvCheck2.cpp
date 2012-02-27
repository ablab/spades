/***************************************************************************
 * Title:          InvCheck2.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <vector>
#include <map>
#include <string>

// mysql stuff
#include "mysql/mysql.h"

// my stuff
#include "utils.h"

class InvCollection {
  std::map<std::string, std::vector> invPos;
};

typedef std::map<std::string, InvCollection> InvMap;

void InitEnv(int argc, char* argv[],
	     std::string &invFile,
	     std::string &alignDBName,
	     ssize_t &diffThreshold);

int main(int argc, char* argv[]) {


  std::string invFileName;
  std::string alignDBName;
  
  ssize_t diffThreshold = 400;

  InitEnv(argc, argv, invFileName, alignDBName, diffThreshold);


}

ssize_t IsTitle(char* line) {
  std::string str(line);
  if (str.find(".fa") >= 0)
    return 1;
  else
    return 0;
}

void ParseInvFile(std::string invFileName, InvMap& invMap) {
  std::ifstream in;
  openck(infFileName, in);
  char buffer[1024]; // why is a 1k block good?
  while(in) {
    in.getline(buffer, 1024);
  }


}


void InitEnv(int argc, char* argv[],
	     std::string &invFile,
	     std::string &alignDBName,
	     std::string &diffThreshold) {

  ssize_t copt;
  ssize_t i;
  while ( (copt=getopt(argc, argv, "t:d:")) != EOF){
    switch(copt) {
    case 't':
      diffThreshold = atoi(optarg);
      continue;
    case 'd':
      alignDBName = optarg;
      continue;
    }
  }
  i = optind;
  if (argc - i < 1) {
    PrintUsage();
  }
  invFile = argv[i];
}

