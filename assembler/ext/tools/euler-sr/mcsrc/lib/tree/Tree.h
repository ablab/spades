/***************************************************************************
 * Title:          Tree.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <sstream>
#include <vector>

template <typename t_data>
class TreeNode {
public:
  virtual TreeNode* Allocate() {
    return new TreeNode<t_data>;
  }
  //  std::set<std::string> leaves;  
  friend std::istringstream & operator>>(std::istringstream &dataText, TreeNode<t_data> *node) {
    dataText >> node->data;
    return dataText;
  }

  void AddChild(TreeNode<t_data>* child) {
    child->parent = this;
    children.push_back(child);
  }
  bool operator==(const TreeNode<t_data> &node) const {
    return this->data == node.data;
  }
  ssize_t size() {
    return children.size();
  }
  std::vector<TreeNode<t_data> *> children;
  TreeNode<t_data> *parent;
  t_data data;
};

template <typename t_node>
class Tree {
public:
  t_node *root;

};

	
#endif
