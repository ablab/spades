/***************************************************************************
 * Title:          CompareAssemblies.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <sstream>

#include "MUMCluster.h"
#include "MUMClusterParser.h"
#include "Contig.h"
#include "ContigGraph.h"

void SummarizeOverlaps(std::vector<Contig> &refContigs,
		       std::vector<Contig> &qryContigs,
		       ContigNameToIndex &refNameToIndex,
		       ContigNameToIndex &qryNameToIndex,
		       std::ofstream &out);

void PrintContigMapGraphs(std::vector<Contig> &refContigs,
			  std::vector<Contig> &qryContigs,
			  ContigNameToIndex &refNameToIndex,
			  ContigNameToIndex &qryNameToIndex,
			  std::string graphName, ssize_t minTreeSize);

void CreateContigMapGraph(std::vector<Contig> &refContigs,
			  std::vector<Contig> &qryContigs,
			  ContigNameToIndex &refNameToIndex,
			  ContigNameToIndex &qryNameToIndex,
			  std::ofstream &out, ssize_t minTreeSize = 0);


int main(int argc, char* argv[]) {
  std::string inFileName, graphName;
  ssize_t minTreeSize;
  if (argc < 3) {
    std::cout << "usage: cmpasm in.cluster graphName" << std::endl;
    std::cout << "just reads a cluster file, nothing else " << std::endl;
    exit(0);
  }
  inFileName = argv[1];
  graphName  = argv[2];
  minTreeSize = 0;
  if (argc > 3) {
    minTreeSize = atoi(argv[3]);
  }
  std::ifstream in;
  in.open(inFileName.c_str());
  MUMClusterFile clusterFile;
  MUMClusterParser::ParseMUMClusterFile(in, clusterFile);

  std::cout << "read " << clusterFile.size() << " clusters " << std::endl;
  ssize_t i;
	//UNUSED// int j;
  ContigNameToIndex refNameToIndex, qryNameToIndex;
  
  std::vector<Contig> refContigs;
  std::vector<Contig> qryContigs;
  CollateClusters(clusterFile, refContigs, qryContigs, refNameToIndex, qryNameToIndex);
  std::cout << "got " << refContigs.size() << " ref contigs and " << qryContigs.size() << " qry contigs " 
	    << std::endl;

  for (i = 0; i < refContigs.size(); i++ ) {
    refContigs[i].FindIntervals();
    refContigs[i].CalculateCoverage();
  }

  std::ofstream graphOut;
  PrintContigMapGraphs(refContigs, qryContigs, refNameToIndex, qryNameToIndex, graphName, minTreeSize);


  std::ofstream summaryOut;
  summaryOut.open("summary");
  SummarizeOverlaps(refContigs, qryContigs, refNameToIndex, qryNameToIndex, summaryOut);
  summaryOut.close();
}


void PrintContigMapGraphs(std::vector<Contig> &refContigs,
			  std::vector<Contig> &qryContigs,
			  ContigNameToIndex &refNameToIndex,
			  ContigNameToIndex &qryNameToIndex,
			  std::string graphName, ssize_t minTreeSize) {
  ForestSet forest;
  BuildForest(refContigs,qryContigs, qryNameToIndex,forest);
  std::ofstream outFile;
  std::stringstream fileNameStrm;
  std::string clusterGraphName;
  ssize_t t,r,q;
  ssize_t treeIndex, qryIndex, qryTreeIndex, qryStartIndex;
  qryStartIndex = refContigs.size();
  for (t = 0; t < forest.size(); t++) {
    // reset the graph name
    if (forest[t].size() > minTreeSize) {
      std::cout << "printing cluster: " << t << " :";
      std::set<ssize_t>::iterator sit, send;
      for (sit = forest[t].begin(); sit != forest[t].end(); sit++) {
	std::cout << " " << *sit;
      }
      std::cout << std::endl;

      
      clusterGraphName = "";
      fileNameStrm.str(clusterGraphName);
      fileNameStrm << graphName << "." << t << ".dot";
      outFile.open(fileNameStrm.str().c_str());
      std::cout << "printing graph " << fileNameStrm.str() << std::endl;
      outFile << "digraph G { " << std::endl;      
      for (r = 0; r < refContigs.size(); r++ ) {
	treeIndex = FindTreeIndex(forest, r);
	//	std::cout << "ref: " << r << " treeIndex: " << treeIndex << " t: " << t << std::endl;
	if (treeIndex == t) {
	  for (q = 0; q < refContigs[r].alignedClusters.size(); q++) {
	    qryIndex = qryNameToIndex[refContigs[r].alignedClusters[q]->qryName] + qryStartIndex;
	    qryTreeIndex = FindTreeIndex(forest, qryIndex);
	    std::cout << "qi: " << qryIndex << " qti: " << qryTreeIndex << std::endl;
	    if (qryTreeIndex == t) {
	      std::stringstream coverageSS;
	      coverageSS << double(refContigs[r].alignedClusters[q]->CalculateCoverage()) / refContigs[r].length;
	      outFile << r << " -> " << qryIndex  
		      << " [label=\"" << coverageSS.str() << "\"];" << std::endl;
	    }
	  }
	}
      }
      for (r = 0; r < refContigs.size(); r++ ) {
	treeIndex = FindTreeIndex(forest, r);
	if (treeIndex == t) 
	  outFile << r << "[ label=\""<< refContigs[r].name 
		  << "," << refContigs[r].length << "\"]" << std::endl;
      }
      for (q = 0; q < qryContigs.size(); q++ ) {
	qryTreeIndex = FindTreeIndex(forest, q + qryStartIndex);
	if (qryTreeIndex == t) 
	  outFile << q + qryStartIndex << "[ label=\""<< qryContigs[q].name 
		  << "," << qryContigs[q].length << "\"]" << std::endl;
      }
    
      outFile << " }" << std::endl;
      outFile.close();
      std::cout<< std::endl;
    }
  }
}

void CreateContigMapGraph(std::vector<Contig> &refContigs,
			  std::vector<Contig> &qryContigs,
			  ContigNameToIndex &refNameToIndex,
			  ContigNameToIndex &qryNameToIndex,
			  std::ofstream &out, ssize_t minTreeSize) {

  out << "digraph G { " << std::endl;
  ssize_t r, q;
  ssize_t qryStartIndex = refContigs.size();
  ssize_t qryIndex;
  ForestSet forest;
  if (minTreeSize > 0) {
    BuildForest(refContigs,qryContigs, qryNameToIndex,forest);
  }

  for (r = 0; r < refContigs.size(); r++ ) {
    if (ValidTreeSize(forest, r, minTreeSize) ) {
      for (q = 0; q < refContigs[r].alignedClusters.size(); q++) {
	qryIndex = qryNameToIndex[refContigs[r].alignedClusters[q]->qryName];
	std::stringstream coverageSS;
	coverageSS << double(refContigs[r].alignedClusters[q]->CalculateCoverage()) / refContigs[r].length;
	out << r << " -> " << qryIndex + qryStartIndex 
	    << " [label=\"" << coverageSS.str() << "\"];" << std::endl;
      }
    }
  }
  
  for (r = 0; r < refContigs.size(); r++ ) {
    if (ValidTreeSize(forest, r, minTreeSize))
	out << r << "[ label=\""<< refContigs[r].name 
	    << "," << refContigs[r].length << "\"]" << std::endl;
  }
  for (q = 0; q < qryContigs.size(); q++ ) {
    if (ValidTreeSize(forest, q + qryStartIndex, minTreeSize))
      out << q + qryStartIndex << "[ label=\""<< qryContigs[q].name 
	  << "," << qryContigs[q].length << "\"]" << std::endl;
  }
  out << " }" << std::endl;
}


void SummarizeOverlaps(std::vector<Contig> &refContigs,
		       std::vector<Contig> &qryContigs,
		       ContigNameToIndex &refNameToIndex,
		       ContigNameToIndex &qryNameToIndex,
		       std::ofstream &out) {
  ssize_t r;
	//UNUSED//  int q;
	//UNUSED//  int qryStartIndex = refContigs.size();
  //UNUSED// ssize_t qryIndex;
  for (r = 0; r < refContigs.size(); r++ ) {
    out << refContigs[r].length << " " << refContigs[r].alignedClusters.size() << std::endl;
  }
}
