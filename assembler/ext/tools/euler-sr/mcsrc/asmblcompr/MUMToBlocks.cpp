/***************************************************************************
 * Title:          MUMToBlocks.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  10/24/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "blocks/BlockGraph.h"

#include "MUMCluster.h"
#include "MUMClusterParser.h"
#include "Contig.h"
#include "ContigGraph.h"


int main(int argc, char* argv[]) {
  std::string inFileName;
  std::string graphFileName;
  if (argc < 2) {
    std::cout << "usage: mum2blocks infile graphfile " << std::endl;
    exit(0);
  }
  inFileName = argv[1];
  graphFileName = argv[2];
  std::ifstream in;
  in.open(inFileName.c_str());
  MUMClusterFile clusterFile;
  MUMClusterParser::ParseMUMClusterFile(in, clusterFile);
  BlockGraph blockGraph;

    ContigNameToIndex refNameToIndex, qryNameToIndex;
  ssize_t refContigIndex, qryContigIndex;
  refContigIndex = -1;
  qryContigIndex = -1;
  
  std::vector<Contig> refContigs;
  std::vector<Contig> qryContigs;
  CollateClusters(clusterFile, refContigs, qryContigs, refNameToIndex, qryNameToIndex);

  ssize_t r, q;
  for (r = 0; r < refContigs.size(); r++) {
    blockGraph.AddVertex(0,refContigs[r].length);
  }

  for (q = 0; q < qryContigs.size(); q++) {
    blockGraph.AddVertex(0,qryContigs[q].length);
  }
  
  ssize_t c, b, p; // cluster, block
  // Now glue together sequences in the graph
  for (r = 0; r < refContigs.size(); r++ ){
    for (c = 0; c < refContigs[r].alignedClusters.size(); c++) {
      MUMCluster *cluster = refContigs[r].alignedClusters[c];
      q = qryNameToIndex[cluster->qryName] + refContigs.size();
      for (b = 0; b < cluster->size(); b++ ) {
	MUMAlignBlock *ab = cluster->blocks[b];
	if (ab->qryStrand == 1) {
	  for (p = 0; p < ab->size(); p++) {
	    blockGraph.GlueSequences(r, 
				     ab->refPos[p] - 1, 
				     ab->refPos[p] + ab->length[p] - 1 - 1,
				     q, 
				     ab->qryPos[p] - 1, 
				     ab->qryPos[p] + ab->length[p] - 1 - 1);
	  }
	}
      }
    }
  }
}


