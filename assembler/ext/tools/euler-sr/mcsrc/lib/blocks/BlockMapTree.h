/***************************************************************************
 * Title:          BlockMapTree.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_MAP_TREE
#define _BLOCK_MAP_TREE

#include "rbtree.h"

#include <ostream>

class Block;

class BlockRef {
  /* 
      BlockRef is used to index into the block graph for 
     easy lookup of vertices.  Since each block contains a list of 
     coordinates in species, the blockref has a reference to the 
     block, and the coordinates of the sequence in the block. 
     This is a bit redundant, but I'll change this when I come up
     with something better.
  */
     
public:
  ssize_t start, end; 
  Block *block;
  BlockRef() { };
  BlockRef(ssize_t pstart, ssize_t pend, Block *pblock) {
    start = pstart;
    end   = pend;
    block = pblock;
  }
  BlockRef &operator=(const BlockRef &bref) { 
		if (this != &bref) {
			start = bref.start;
			end   = bref.end;
			block = bref.block;
		}
    return *this;
  }
  bool operator==(const BlockRef &b) const {
    // Not a transitive operation, tests 
    // for b included in this block
    return (start <= b.start &&
	    end >= b.end);
  }
  
  bool operator>(const BlockRef &b) const {
    return (start >= b.start);
  }
  
  bool operator<(const BlockRef &b) const {
    return (start < b.start);
  }
};


typedef RBTreeNode<BlockRef> BlockNode;

class BlockMapTree : public RBTree<BlockRef> {
  BlockNode* FindEnclosingBlockRef(ssize_t pos, BlockNode *node=NULL);
public:
  void Print(std::ostream &out, BlockNode *node=NULL, std::string padding="");
  bool RemoveBlockRef(ssize_t pos);
  Block* FindEnclosingBlock(ssize_t pos);
  BlockNode *begin() { return TreeMinimum(); }
};

class BlockMapIterator {
public:
  BlockNode *cur;
  BlockMapTree *tree;
  BlockMapIterator() { cur = NULL; tree = NULL; }
  BlockMapIterator &operator=(const BlockMapIterator &rhs) {
		if (this != &rhs) {
			cur = rhs.cur;
			tree = rhs.tree;
		}
		return *this;
  }

  Block* Data() {
    if (cur != NULL)
      return cur->data.block;
		return (Block *) NULL;
  }
  
  BlockMapIterator & Next() {
    cur = tree->Successor(cur);
    return *this;
  }
  
  BlockMapIterator & operator++() {
    Next();
    return *this;
  }

  ssize_t NotNil() {
    return (cur != tree->nil && cur != NULL);
  }
};
// Define the datatype for mapping from block position to edge in graph


#endif
