/***************************************************************************
 * Title:          BlockMap.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_MAP_
#define _BLOCK_MAP_

#include "rbtree.h"
#include "BlockMapTree.h"
#include "Block.h"

#include <ostream>
#include <iostream>


class BlockMap {
  BlockMapTree  blockTree;
  BlockNode* FindNode(ssize_t pos, BlockNode *node);
public:
  Block* begin() {
    BlockNode *bn;
    bn = blockTree.TreeMinimum();
    if (bn != NULL)
      return bn->data.block;
  };
  Block* FindEnclosingBlock(ssize_t pos);
  void AddReference(Block *block, ssize_t index);

  // Remove the reference that contains this pos
  void RemoveReference(ssize_t pos);
  bool RemoveReference(Block *block, ssize_t index);

  void BeginIterator(BlockMapIterator &bit) {
    bit.cur  = blockTree.begin();
    bit.tree = &blockTree; // expose some data, oh well.
  }
  void Print(std::ostream &out) {blockTree.Print(out);}
  
  ssize_t size() { return blockTree.size(); }
};

#endif
