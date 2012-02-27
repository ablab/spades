/***************************************************************************
 * Title:          BlockGraph.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_GRAPH_H
#define _BLOCK_GRAPH_H

#include "Block.h"
#include "BlockMap.h"
#include "BlockPath.h"

#include <vector>

class BlockGraph {
  // This class implements the data structure to hold a representation of a 
  // set of alignments of genomes as a block graph.  
  // Each edge is a sequence that maps to one or more genomes.  The 
  // length of each edge is the number of nucleotides that it maps to 
  // in one genome.  The edge is a gap-free alignment
  // (the edge maps to the same number of nucleotides in every genome). 

public:
  // Index into the graph via a map that hashes according to the starting position
  // of the edge in the reference genome.

  std::vector<BlockMap*> vertexMaps;
  BlockGraph() {
    // do nothing for now in the constructor.
  }

  // Create an unconnected sequence from the graph
  ssize_t AddVertex(ssize_t start, ssize_t end);

  // Join the paths corresponding to two sequences into one.
  void GlueSequences(ssize_t refStrand, ssize_t startRef, ssize_t endRef,
		     ssize_t qryStrand, ssize_t startQry, ssize_t endQry);

  // Find the vertex following this one.  I might just use the 
  // rb tree for this, but for now am using the graph.
  Block* FindNext(Block *block, ssize_t seq, ssize_t pos, Block *&next, ssize_t &s);

  // Given two paths in this graph, merge them into one.  The coordinates
  // of each path are taken to be the coordinates of the sequence in
  // the path.
  void MergePaths(BlockPath &refPath, BlockPath &qryPath);
  // Locate a path in the graph
  void FindPath(ssize_t seq, ssize_t start, ssize_t end, BlockPath &path);

  // Take a vertex and split it into two, or three.
  void SplitVertex(ssize_t sequence, ssize_t start, ssize_t end);
  void SplitVertex(Block *&vertex,  ssize_t pos, ssize_t seq, Block *&in, Block *&out);
  void SplitVertex(Block *&vertex, ssize_t seq, ssize_t start, ssize_t end);

  void PrintPath(ssize_t map, ssize_t start, ssize_t end, std::ostream &out);

  void ResetVisited();
  void PrintGraph(std::ostream &out);
};


#endif
