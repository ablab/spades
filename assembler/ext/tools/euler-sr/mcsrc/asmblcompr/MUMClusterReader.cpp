/***************************************************************************
 * Title:          MUMClusterReader.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "MUMCluster.h"
#include "MUMClusterParser.h"
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <map>

int main(int argc, char* argv[]) {
  std::string inFileName;
  if (argc < 2) {
    std::cout <<  "usage: mumcr in.cluster" << std::endl;
    std::cout << "just reads a cluster file, nothing else " << std::endl;
    exit(0);
  }
  inFileName = argv[1];
  std::ifstream in;
  in.open(inFileName.c_str());
  MUMClusterFile clusterFile;
  MUMClusterParser::ParseMUMClusterFile(in, clusterFile);

  std::cout << "read " << clusterFile.size() << " clusters " << std::endl;
  ssize_t i;
	//UNUSED// int j;
  for (i = 0; i < clusterFile.size(); i++) {
    std::cout << "cluster: " << i << " size: " << clusterFile.clusters[i]->size() << std::endl;
  }
}
