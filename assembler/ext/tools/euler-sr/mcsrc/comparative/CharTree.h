/***************************************************************************
 * Title:          CharTree.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CHARTREE_H_
#define CHARTREE_H_

#include <string>
#include <set>
#include <vector>
#include <ostream>

#include "tree/Tree.h"


class CharNode : public TreeNode<std::string> {
public:
  ssize_t treeNumber;
  ssize_t value;
  std::set<std::string> leaves;
  std::vector<std::string> locations;
  CharNode() {
    treeNumber = -1;
    value = -1;
  }
  virtual CharNode* Allocate() {
    return new CharNode;
  }
  void Print(ssize_t &counter, std::ostream& out);
  void NewickPrint(std::ostream& out);
};

class CharTree : public Tree<CharNode> {
  void StoreInternalList(CharNode *node, ssize_t &number);
public:
  void StoreInternalList() {
    ssize_t number = 0;
    StoreInternalList(static_cast<CharNode*>(root), number);
  }
  ssize_t GetListSimilarity(CharNode *node, std::set<std::string> &qry);
  void LocateClade(std::set<std::string> &cladeSet, 
		   CharNode *node,
		   CharNode *&minCladeNode,
		   ssize_t &minScore);

  void Print(std::ostream &out);
  void NewickPrint(std::ostream &out);

  CharNode* FindNode(CharNode *query, CharNode *node= NULL) {
    if (node == NULL)
      node = this->root;
    // DFS search for child;
    ssize_t i;
    if (*query == *node)
      return node;
    
    CharNode *result;
    for (i = 0; i < node->children.size(); i++) {
      if ((result = FindNode(query, static_cast<CharNode*>(node->children[i]))) != NULL) {
	return result;
      }
    }
    return NULL;
  }

};

#endif
