/***************************************************************************
 * Title:          CharTree.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <fstream>
#include <list>

#include "CharTree.h"
#include "tree/NewettTree.h"
#include "utils.h"


void CharNode::NewickPrint(std::ostream &out) {
  ssize_t i;
  out << "(";
  if (children.size() > 0) {
    for (i = 0; i < children.size()-1; i++) {
      static_cast<CharNode*>(children[i])->NewickPrint(out);
      out << ",";
    }
    static_cast<CharNode*>(children[i])->NewickPrint(out);
  }
  else 
    out << data;
  out << ")";
}

void CharNode::Print(ssize_t &counter, std::ostream &out) {

  ssize_t i;
  if ( treeNumber == -1 ) 
    treeNumber = ++counter;
    
  for (i = 0; i < children.size(); i++) {
    if (static_cast<CharNode*>(children[i])->treeNumber == -1)
      static_cast<CharNode*>(children[i])->treeNumber = ++counter;
    out << treeNumber << " -> " << static_cast<CharNode*>(children[i])->treeNumber 
	<< "; " << std::endl;
  }
  out << treeNumber << "[ shape= box ];" << std::endl;
  if (value != -1)
    out << treeNumber << "[ label= \"" << data << " " << value << "\\n ";
  else
    out << treeNumber << "[ label= \"" << data << "\\n ";
  if (locations.size() > 0)
    out << "(";
  for (i = 0; i < ((ssize_t)locations.size())-1; i++) 
    out << " " << locations[i] << ",";
  if (locations.size() > 0) 
    out << " " << locations[i] << ")";
    
  out << "\"]; " << std::endl;
  ++counter;
  for (i = 0; i < children.size(); i++) {
    static_cast<CharNode*>(children[i])->Print(counter, out);
  }
}

void CharTree::NewickPrint(std::ostream &out) {
  static_cast<CharNode*>(root)->NewickPrint(out);
  out << ";" << std::endl;
}

void CharTree::Print(std::ostream &out) {
  if (root == NULL) 
    return;
  ssize_t counter = 0;
  out << "digraph G { " << std::endl;
  out << "size=\"8,10\";" << std::endl;
  static_cast<CharNode*>(root)->treeNumber = counter;
  static_cast<CharNode*>(root)->Print(counter, out);
  out << " }" << std::endl;
}

void CharTree::StoreInternalList(CharNode *node, ssize_t &number) {
  /* 
     This builds a list of all the data at the leaves.  
     On a big tree it's approximately logn*n additional size (if it's balanced
     binary).
  */
  ssize_t i, j;
  if (node->size() == 0) 
    node->leaves.insert(node->data);

  node->treeNumber = number;

  for (i = 0; i < node->size(); i++) {
    StoreInternalList(static_cast <CharNode*>(node->children[i]), ++number);
    (static_cast<CharNode*>(node))->
      leaves.insert((static_cast<CharNode*>(node->children[i]))->leaves.begin(),
		    (static_cast<CharNode*>(node->children[i]))->leaves.end());
  }
}


ssize_t CharTree::GetListSimilarity(CharNode *node, std::set<std::string> &qry) {

  ssize_t qryNotTarget, targetNotQry;
  std::list<std::string> symmDiff, excludedDiff;
  std::front_insert_iterator<std::list<std::string> > diff(symmDiff); 
  std::front_insert_iterator<std::list<std::string> > excluded(excludedDiff);  
  
  std::set<std::string>::iterator nodeit, qryit;
  /*
  std::cout << "node:";
  for (nodeit = node->leaves.begin(); nodeit != node->leaves.end(); ++nodeit) 
    std::cout << " " << *nodeit ;
  std::cout << std::endl;
    
  std::cout << "qry:";
  for (qryit = qry.begin(); qryit != qry.end(); ++qryit) 
    std::cout << " " << *qryit;
  std::cout << std::endl;
  */
  std::set_difference(qry.begin(), qry.end(), 
		      node->leaves.begin(), node->leaves.end(),
		      excluded);

  if (excludedDiff.size() > 0) {
    std::list<std::string>::iterator lit;
    //    std::cout << " excluded " << std::endl;
    return 100000000; // This node does not contain everything.
  }
  else {
    //    std::cout << "node contains query " << std::endl;
    // This node contains all of the query set, try to find how similar 
    // this node is to the query set.
    std::set_symmetric_difference(node->leaves.begin(), node->leaves.end(),
				  qry.begin(), qry.end(),
				  diff);
    /*
    std::cout << "node: " << node->treeNumber<< " " 
	      << node->leaves.size() << " qry: " << qry.size() << " diff: " << symmDiff.size() << std::endl;
    */
    return symmDiff.size();
  }
}


void CharTree::LocateClade(std::set<std::string> &cladeSet, 
			   CharNode *node,
			   CharNode *&minCladeNode,
			   ssize_t &minScore) {
  
  ssize_t score;
  score = GetListSimilarity(node, cladeSet);
  if (score < minScore) {
    minScore = score;
    minCladeNode = node;
  }
  ssize_t i;
  for (i = 0; i < node->size(); i++) {
    LocateClade(cladeSet, static_cast<CharNode*>(node->children[i]),
		minCladeNode, minScore);
  }
}
