/***************************************************************************
 * Title:          CompareContigs.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <map>

#include "MUMCluster.h"
#include "MUMClusterParser.h"
#include "Contig.h"
#include "ContigGraph.h"


typedef std::map<std::string, ssize_t> ContigNameToIndex;


void SummarizeOverlaps(std::vector<Contig> &refContigs,
		       std::vector<Contig> &qryContigs,
		       ContigNameToIndex &refNameToIndex,
		       ContigNameToIndex &qryNameToIndex,
		       std::ofstream &out);

void CompareOneToOne(Contig &refContig, Contig &qryContig, double &refCoverage, double &qryCoverage) {
  refCoverage = ((double) refContig.covered )/ refContig.length;
  qryCoverage = ((double) qryContig.covered )/ qryContig.length;
}

int main(int argc, char* argv[]) {
  std::string inFileName, graphName;
  if (argc < 2) {
    std::cout <<  "usage: mumcr in.cluster " << std::endl;
    std::cout << "just reads a cluster file, nothing else " << std::endl;
    std::cout << "will print some interesting statistics on an assembly " << std::endl;
    exit(0);
  }
  inFileName = argv[1];
  std::ifstream in;
  in.open(inFileName.c_str());
  MUMClusterFile clusterFile;
  MUMClusterParser::ParseMUMClusterFile(in, clusterFile);

  std::cout << "read " << clusterFile.size() << " clusters " << std::endl;
  //UNUSED// int i, j;
  ContigNameToIndex refNameToIndex, qryNameToIndex;
  ssize_t refContigIndex, qryContigIndex;
  refContigIndex = -1;
  qryContigIndex = -1;
  
  std::vector<Contig> refContigs;
  std::vector<Contig> qryContigs;
  Contig emptyContig;

  CollateClusters(clusterFile, refContigs, qryContigs, refNameToIndex, qryNameToIndex);

  std::cout << "got " << refContigs.size() << " ref contigs and " 
	    << qryContigs.size() << " qry contigs " << std::endl;


  ForestSet forest;  // it's just a bag of nodes, with no edges,but
		     // we'll call it a forest anyway  
  
  BuildForest(refContigs,qryContigs, qryNameToIndex,forest);

  // Summarize the resulting graphs.
  
  ssize_t n1to1 = 0;
  ssize_t t;
  std::vector<ssize_t> refIndices, qryIndices;
  double refCoverage, qryCoverage;
  for (t = 0; t < forest.size(); t++) {
    if (forest[t].size() == 2) {
      refIndices.clear();
      qryIndices.clear();
      ChopTree(forest[t], refContigs.size(), refIndices, qryIndices);
      assert(refIndices.size() == 1);
      assert(qryIndices.size() == 1);
      CompareOneToOne(refContigs[refIndices[0]], qryContigs[qryIndices[0]], refCoverage, qryCoverage);
      std::cout << refIndices[0] << ": " << refCoverage << " " << qryIndices[0] << ": " << qryCoverage<<std::endl;
      ++n1to1;
    }
  }
  std::cout << "found: " << n1to1 << " 1-1 contigs " << std::endl;

  std::ofstream summaryOut;
  summaryOut.open("summary");
  SummarizeOverlaps(refContigs, qryContigs, refNameToIndex, qryNameToIndex, summaryOut);
  summaryOut.close();
}


		     



void SummarizeOverlaps(std::vector<Contig> &refContigs,
		       std::vector<Contig> &qryContigs,
		       ContigNameToIndex &refNameToIndex,
		       ContigNameToIndex &qryNameToIndex,
		       std::ofstream &out) {
  ssize_t r;
	//UNUSED// int q;
  //UNUSED// int qryStartIndex = refContigs.size();
  //UNUSED// int qryIndex;
  for (r = 0; r < refContigs.size(); r++ ) {
    out << refContigs[r].length << " " << refContigs[r].alignedClusters.size() << std::endl;
  }
}
