/***************************************************************************
 * Title:          MUMClusterParser.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <assert.h>
#include <stdlib.h>
#include "MUMClusterParser.h"

ssize_t MUMClusterParser::ParseMUMClusterFile(std::ifstream &inFile,
					  MUMClusterFile &clusterFile) {
  inFile >> clusterFile.refFileName >> clusterFile.qryFileName;
  inFile >> clusterFile.method;
  //  std::cout << "got method: " << clusterFile.method << std::endl;
  inFile.get();
  MUMCluster *cluster;
  while (ParseCluster(inFile, cluster) and cluster != NULL) {
    clusterFile.clusters.push_back(cluster);
  }
  if (cluster == NULL)
    return 0;
  else {
    clusterFile.clusters.push_back(cluster);
    return cluster->size();
  }
}
					  

ssize_t MUMClusterParser::ParseCluster(std::ifstream &inFile, MUMCluster *&cluster) {
  cluster = NULL;
  if (inFile.good() and inFile.peek() == '>') {
    cluster = new MUMCluster;
		//    char c = inFile.get(); // discard the '>'
		inFile.get(); // discard the '>'
    //    std::cout << "got char: " << c << std::endl;
    inFile >> cluster->refName >> cluster->qryName >> cluster->refLen >> cluster->qryLen;
    inFile.get();
    ParseAlignBlocks(inFile, cluster);
  }
  return inFile.good();
}


ssize_t MUMClusterParser::ParseAlignBlocks(std::ifstream &inFile, 
				       MUMCluster *cluster) {
  ssize_t refPos, qryPos, length;
  std::string refGap, qryGap;
  // Initialize the block to assume we haven't read it.
  assert(cluster != NULL);

  MUMAlignBlock *block;
  block = NULL;

  // If a align block exists, read it 
  while ( inFile and 
	  inFile.peek() != '>' and 
	  inFile.peek() != '\0') {
    inFile >> refPos >> qryPos;
    if (inFile.peek() == '\n') {
      // The two entries that mark the start of a new block
      block = new MUMAlignBlock;
      cluster->blocks.push_back(block);
      block->refStrand = refPos;
      block->qryStrand = qryPos;
    }
    else {
      // Read in an alignment of the block
      inFile >> length >> refGap >> qryGap;
      block->refPos.push_back(refPos);
      block->qryPos.push_back(qryPos);
      block->length.push_back(length);
      if (refGap == "-")
	block->refGap.push_back(-1);
      else
	block->refGap.push_back(atoi(refGap.c_str()));
      
      if (qryGap == "-")
	block->qryGap.push_back(-1);
      else
	block->qryGap.push_back(atoi(qryGap.c_str()));
    }
    inFile.get();
  }

  char pstatus = inFile.peek();
  bool fstatus = inFile.eof();
  return !fstatus  and pstatus != '>';
}


