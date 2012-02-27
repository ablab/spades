/***************************************************************************
 * Title:          BlockMapTree.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BlockMapTree.h"

#include "Block.h"

void BlockMapTree::Print(std::ostream &out, 
			 BlockNode *node, 
			 std::string padding) {
  // don't even try to print if there is nothing here
  if (root == NULL || root == nil) return; 

  // default parameters starts on root
  if (node == NULL) 
    node = root;
  
  out << padding << *(node->data.block) << "{" << node << "} (" << node->leftChild << ", " 
      << node->rightChild << ")" << std::endl;

  if (node->leftChild != nil)   // only recurse on non-null nodes
    Print(out, node->leftChild, padding+" l ");
  else
    out << padding + " " << "nil" << std::endl;
  if (node->rightChild != nil)
    Print(out, node->rightChild, padding + " r ");
  else
    out << padding + " " << "nil" << std::endl;
}
 
BlockNode* BlockMapTree::FindEnclosingBlockRef(ssize_t pos, BlockNode *node) {
  if (node == NULL)
    node = root;
  if (node == NULL) 
    return NULL;

  if (node->data.start <= pos &&
      node->data.end >= pos)
    return node;
  else {
    if (node->data.start > pos) {
      if (node->leftChild == NULL || node->leftChild == nil)
				return NULL;
      else
				return FindEnclosingBlockRef(pos, node->leftChild);
    }
    if (node->data.end < pos) {
      if (node->rightChild == NULL || node->rightChild == nil)
				return NULL;
      else
				return FindEnclosingBlockRef(pos, node->rightChild);
    }
  }
	assert(0);
	return node; // Never reaches here.  Quiet compiler warnings.
}

Block* BlockMapTree::FindEnclosingBlock(ssize_t pos) {
  BlockNode* blockNode;
  blockNode = FindEnclosingBlockRef(pos);
  if (blockNode != NULL) 
    return blockNode->data.block;
  else
    return NULL;
}

bool BlockMapTree::RemoveBlockRef(ssize_t pos) {
  BlockNode* blockNode;
  blockNode = FindEnclosingBlockRef(pos);
  if (blockNode != NULL) {
    Delete(blockNode);
    delete blockNode;
    return true;
  }
  else {
    return false;
  }
}
