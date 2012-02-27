/***************************************************************************
 * Title:          BlockMap.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <assert.h>
#include "BlockMap.h"

Block* BlockMap::FindEnclosingBlock(ssize_t pos) {
  return blockTree.FindEnclosingBlock(pos);
}

void BlockMap::AddReference(Block *block, ssize_t index) {
  // parameters: block, make the map reference this block.
  //             index: 
  ssize_t bIndex;
  bIndex = block->GetSequenceIndex(index);
  assert(block->size() > bIndex);
  std::cout << "adding reference to " << index << "(" << bIndex << ")"
	    << " start: " << block->start[bIndex] 
	    << " end: " << block->end[bIndex] << std::endl;


  BlockNode *node = blockTree.Insert(BlockRef(block->start[bIndex], 
					      block->end[bIndex], block));

  //  blockTree.Print(std::cout);
  block->ref[index] = node;
  std::cout << "add index " << bIndex << " to " << block << std::endl;
}

bool BlockMap::RemoveReference(Block* block, ssize_t index) {
  // new way of removing nodes
  //  BlockNode *newNode; 
  //  newNode = blockTree.Find(BlockRef(node->data.block->start[0], node->data.block->end[0], node->data.block));
  bool result;
  BlockNode *node;
  node = block->ref[index];
  if (node != NULL) {
    std::cout << "deleting index " << index << " from " << block << std::endl;

  /*
    std::cout << "tree before " << std::endl;
    blockTree.Print(std::cout);
    std::cout << "deleting tree node " << node << std::endl;
  */
    result = blockTree.Delete(node);
    block->ref[index] = NULL;
    /*
      std::cout << "tree afterwards " << std::endl;
      blockTree.Print(std::cout);
    */
    delete node;
    return result;
  }
  else 
    return 1;
}

void BlockMap::RemoveReference(ssize_t pos) {
  blockTree.RemoveBlockRef(pos);
}

