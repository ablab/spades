/***************************************************************************
 * Title:          Block.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/28/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_H_
#define _BLOCK_H_

#include <vector>
#include <map>
#include <stdio.h>
#include <iostream>
#include <ostream>
#include <assert.h>
#include "BlockMapTree.h"

class Block {
  // A block is an ungapped sequence shared between 1 or more sequences.  It contains 
  // the coordinates of the seuqence in each species, and the label of
  // the species that each pair of coordinates corresponds to.  
public:
  // A block represents a contiguous stretch in a collection of
  // sequences indexed 0 .. N.  Store the indices of the sequences
  // that map to this block.

  // use this attribute for bfs type searches through the graph.
  char marked;
  ssize_t number;

	//////// TODO: Here, marked is character and 0=unmarked, 1=marked
	//////// Whereas for GraphVertex and GraphEdge,
	//////// it's a one-bit enum, 0=Marked, 1=NotMarked

  // There is one map per sequence. The map is a collection of
  // blocknodes that each point to a block.  Each block references all
  // block nodes that reference it so that updates may be made to all
  // maps that reference this block if this block is changed.

  std::map<ssize_t, BlockNode*> ref;

  // The indices of the sequences referenced by this block
  std::vector<ssize_t> sequence;

  // Store the coordinates of each block.  I'm not sure yet what
  // convention to use for location information.  I'll probably use
  // first position = 1, and sign indicates strand.  Reverse
  // complement strand is indexed on the coordinate system of the
  // forward strand.
  std::vector<ssize_t> start, end;
  std::vector<Block*> in, out;

  Block() {
    marked = 0;
    number = 0;
  }
  Block(ssize_t startPos, ssize_t endPos, ssize_t seq) {
    std::cout << "creating block for " << seq << " in block " << this << std::endl;
    sequence.push_back(seq);
    start.push_back(startPos);
    end.push_back(endPos);
    ref[seq] = NULL;
    marked = 0;
    number = 0;
  }

  ssize_t size() const {
    assert(sequence.size() == start.size() && 
	   start.size() ==end.size());
    return sequence.size();
  }

  void AddCoords(ssize_t qryStrand, ssize_t qryStart, ssize_t qryEnd);

  // Take a window that corresponds to a vertex, and 
  // split it into some region before the window, 
  // an intersection with the window, and a region after 
  // the window.  It's possible that the region before (in)
  // and the region after (out) are null.
  void Split(ssize_t start, ssize_t end, ssize_t seq,
	     Block *&in, Block *&intersect, Block *&out);

  // Take one vertex and split it into two.
  void Split(ssize_t pos, ssize_t seq, Block *&inVertex, Block *&outVertex);

  // Merge attributes (in / out edges, sequence coordinates)
  // into this vertex.
  void Merge(Block *vertex);

  ssize_t GetSequenceIndex(ssize_t seq, ssize_t startIndex=0);
  ssize_t GetSequenceIndexPos(ssize_t seq, ssize_t pos, ssize_t startIndex=0);
  ssize_t GetRefPos(ssize_t pos, ssize_t seqIndex);
  void Print(std::ofstream &out);
};

std::ostream &operator<<(std::ostream &out, const Block &b);

#endif
