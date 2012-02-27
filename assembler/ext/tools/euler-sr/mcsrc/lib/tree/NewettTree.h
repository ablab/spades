/***************************************************************************
 * Title:          NewettTree.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NEWETT_TREE_H_
#define NEWETT_TREE_H_
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "Tree.h"


void ReadData(std::ifstream &in, std::string &data);
void Advance(std::ifstream &in);  

template <typename T>
void ReadSubTree(std::ifstream &in, TreeNode<T> *localRoot) {
  std::string data;
  std::istringstream datastr;
  char next;
  Advance(in);
  next = in.get();
  //  std::cout << "got: " << next << std::endl;
  if (next != '(' and
      next != ')') {
    std::cout << "invalid tree format near" ;
    std::string errorstr;
    in >> errorstr;
    std::cout << errorstr << std::endl;
    exit(0);
  }
  ReadData(in, data); // advances to next ')'
  datastr.str(data);
  
  // store the contents of this node.
  datastr >> localRoot;
  in >> std::ws;
  next = in.peek();
  while (next != ')' and next != ';' and  next != '\0' and in.good()) {
    if (next == ',') {
      in.get();
      next = in.peek();
    }
    if (next == '(') {
      TreeNode<T>*child = localRoot->Allocate(); //new Node<T>;
      localRoot->AddChild(child);
      ReadSubTree(in, child);
    }
    else {
      if (next != ')') {
	std::cout << " invalid tree format near: '" << next << "' " ;
	std::string errorstr;
	in >> errorstr;
	std::cout << errorstr << std::endl;
	std::cout << "expected ')' " << std::endl;
      }
    }
    next = in.peek();
  }
  in.get();
  Advance(in);
}

template<typename T>
void ReadTree(std::ifstream &in, Tree<T> &tree) {
  if (in.peek() != '(') {
    std::cout << "badly formmatted tree, should start with '(' "
	      << std::endl;
    exit(1);
  }
  tree.root =  new T;
  ReadSubTree(in, tree.root);
}

template<typename T>
void PrintTreeNode(std::ostream &out, TreeNode<T> *node) {
  std::cout << node->data << std::endl;
  ssize_t i;
  std::cout << ">>>>" << std::endl;
  for (i = 0; i < node->size(); i++) {
    PrintTreeNode(out, node->children[i]);
  }
  std::cout << "<<<" << std::endl;
}


template<typename T>
void PrintTree(std::ostream &out, Tree<T> &tree) {
  PrintTreeNode(out, tree.root);
}
#endif
