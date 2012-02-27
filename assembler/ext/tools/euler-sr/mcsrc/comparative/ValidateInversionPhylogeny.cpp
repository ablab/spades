/***************************************************************************
 * Title:          ValidateInversionPhylogeny.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <sstream>

#include <set>

#include "CharTree.h"
#include "tree/NewettTree.h"
#include "utils.h"
#include "InversionChars.h"
#include "InversionBins.h"

typedef std::map<std::string, CharNode *> NodeMap;

ssize_t verbose;
typedef std::map<std::string, ssize_t> CharMap;
ssize_t CheckPhylogeny(CharNode* cNode, 
		   CharMap &charMap, 
		   ssize_t &value, ssize_t &branchFound, ssize_t &violated, 
		   std::string &branchData, std::string &violatingNode);

ssize_t BuildNodeMap(StringVector &species, CharTree &tree, NodeMap &nodeMap);

int main(int argc, char* argv[]) {
  std::string binFileName, treeFileName, graphFileName;
  StringSet ignoreSpecies;
  NodeMap nodeMap;
  std::string rejectedFileName, fixedFileName;
  if (argc < 4) {
    std::cout << "usage: validphylo binFileName treeFileName graphFileName [-v] [-r rejectedfile] [-f fixedFileName]" 
	      << " [ignorespec1] [ignorespec2] ..." << std::endl;
    exit(1);
  }

  binFileName   = argv[1];
  treeFileName  = argv[2];
  graphFileName = argv[3];
  rejectedFileName = "";
  fixedFileName = "";
  int argi = 4;
  verbose = 0;
  while (argi < argc) {
    if (strcmp(argv[argi], "-v") == 0)
      ++verbose;
    else {
      ignoreSpecies.insert(std::string(argv[argi]));
    }
    if (strcmp(argv[argi], "-r") == 0) {
      ++argi;
      assert(argi < argc);
      rejectedFileName = argv[argi];
    }
    if (strcmp(argv[argi], "-f") == 0) {
      ++argi;
      assert(argi < argc);
      fixedFileName = argv[argi];
    }
    ++argi;
  }

  // Read in the strings to validate
  StringVector species;
  std::vector<ssize_t> startPos, endPos;
  BinMap binnedInversions;
  ReadBinFile(binFileName, species, binnedInversions);

  InversionList consensus;
  std::vector<ssize_t> mult;
  GetBinConsensus(binnedInversions, species.size(), startPos, endPos, consensus, mult);
  
  CharMap charValues;

  CharTree cTree;
  std::ifstream treeIn;
  std::set<std::string> setA, setB;
  std::vector<std::string> *invSpecP;
  openck(treeFileName, treeIn);
  ReadTree<CharNode>(treeIn, cTree); 

  // Store a map of the tree edges so they can be accessed easily.
  BuildNodeMap(species, cTree, nodeMap);

  // Store a list of the children of every node.
  // not needed  cTree.StoreInternalList();
  //

  ssize_t value, branchFound, violated;
  std::string branchData;
  // validate all of the validated inversions.
  ssize_t bin, inv;
  InversionList::iterator invIt, prevIt;
  std::string violatingNode;
  for (bin = 0; bin < binnedInversions.size(); bin++) {
    invIt = binnedInversions[bin].begin(); 
    while (binnedInversions[bin].size() > 0 and
	   invIt != binnedInversions[bin].end() ) {
      // Store the value of the chars per species.
      if (ignoreSpecies.find((*invIt)->species) != ignoreSpecies.end())
	continue;

      for (inv = 0; inv < (*invIt)->inversions.size(); ++inv) {
	if (ignoreSpecies.find(species[inv]) == ignoreSpecies.end())
	  charValues[species[inv]] = (*invIt)->inversions[inv];
	else
	  charValues[species[inv]] = 2;
	nodeMap[species[inv]]->value = charValues[species[inv]];
      }

      value = 0;
      branchFound = 0;
      violated = 0;
      branchData = "";
      violatingNode = "";

      CheckPhylogeny(cTree.root, charValues, value, 
		     branchFound, violated, branchData, violatingNode);
      if (violated == 1) {
	std::cout << "checking phylogeny got: " << (*invIt)->startPos 
		  << " " << (*invIt)->endPos << " " 
		  << branchFound << " " << branchData << " " << violated 
		  << " " << violatingNode << std::endl;
	std::ofstream rejectedFile;
	std::stringstream sstr;
	sstr << rejectedFileName << "_" << (*invIt)->startPos << ".dot";
	std::string rejectedTreeFileName = sstr.str();
	openck(rejectedTreeFileName, rejectedFile);
	cTree.Print(rejectedFile);
	rejectedFile.close();
	prevIt = invIt;
	++invIt;
	if (invIt == binnedInversions[bin].end()) {
	  binnedInversions[bin].erase(prevIt);
	  invIt = binnedInversions[bin].end();
	}
	else {
	  binnedInversions[bin].erase(prevIt);
	}
      }	
      else {
	++invIt;
      }
    }
  }
  if (fixedFileName != "") {
    std::ofstream fixedOut;
    openck(fixedFileName, fixedOut);
    PrintBins(species, binnedInversions, fixedOut);
    fixedOut.close();
  }
  return 0;
}

ssize_t BuildNodeMap(StringVector &species, CharTree &tree, NodeMap &nodeMap) {
  ssize_t s;
  CharNode query;
  for (s = 0 ; s < species.size(); s++) {
    query.data = species[s];
    nodeMap[species[s]] = static_cast<CharNode*>(tree.FindNode(&query));
    //    nodeMap[species[s]] = tree.FindNode(dynamic_cast<CharNode*>(&query));
  }
}


ssize_t CheckPhylogeny(CharNode* cNode, 
		   CharMap &charMap, 
		   ssize_t &value, ssize_t &branchFound, ssize_t &violated, 
		   std::string &branchData, std::string &violatingNode) { 
  // If this is s leaf, then there is no phylogeny here, the value
  // of the subtree is the value of the leaf
  if (cNode->children.size() == 0) {
    if (verbose)
      std::cout << "leaf: " << cNode->data << " has value: " << charMap[cNode->data]<< std::endl;
    value = charMap[cNode->data];
    return 0;
  }

  // There were leaves here.  Find their values.

  std::vector<ssize_t> childValues, branches;  // sheesh, they're just getting worse

  ssize_t i;
  ssize_t childValue;
  ssize_t childBranch;
  ssize_t numBranchingChildren;
  ssize_t numBranches;
  numBranches = 0;
  for (i = 0; i < cNode->children.size() and !violated; i++) {
    childBranch = 0;
    numBranches += CheckPhylogeny(static_cast<CharNode*>(cNode->children[i]), 
				  charMap, childValue, childBranch, violated, 
			          branchData, violatingNode);
    childValues.push_back(childValue);
    branches.push_back(childBranch);
    if (verbose)
      std::cout << "internal node " << cNode->data << " has child with value: " 
		<< childValue << " " << childBranch << " " << violated << std::endl;
  }

  // If more than one child causes a branch there is a perfect phylogeny violation
  if (verbose) 
    std::cout << "internal node " << cNode->data << " has " << numBranches << " branching vertices " << std::endl;
  if (numBranches > 1) {
    violatingNode = cNode->data;
    violated = 1;
    return numBranches;
  }

  // If only one of the children was a branch here, that's ok (good), mark that this was the branch
  // point
  if (numBranches == 1) {
    // There was only one branch here, but there may be branches further up 
    // in the tree that cause a violation.  The value of the branching edge is 
    // not known, but the value of the other nonbranching edge 
    value = 2; // start off with value unkown (there may be only one edge here)
    for (i = 0; i < cNode->children.size(); i++) {
      if (branches[i] != 1) {
	value = childValues[i];
	break;
      }
    }
  }


  // Determine if there was a branch at this vertex.
  ssize_t zeroFound = 0;
  ssize_t oneFound  = 0;

  for (i = 0; i < cNode->children.size(); i++) {
    if (childValues[i] == 0)
      zeroFound = 1;
    if (childValues[i] == 1)
      oneFound  = 1;
  }

  if (zeroFound and oneFound) {
    // There was a branch, increase the # branches by 1
    branchFound = 1;
    branchData = cNode->data;
    if (verbose) 
      std::cout << "internal node " << cNode->data << " is brancing " << std::endl;
    // Consider this to be unknown, and keep searching for branches higher up in the tree
    // since they are problematic as well.
    value = 2;
    return numBranches + 1; 
  }
  else {
    // No branch at this vertex.  Record the value as one 
    // of the known values of children.
    branchFound = numBranches;
    if (zeroFound) {
      if (verbose) 
	std::cout << cNode->data << " setting value: 0 " << std::endl;
      value = 0;
    }
    else if (oneFound) {
      if (verbose)
	std::cout << cNode->data << " setting value 1 " << std::endl;
      value = 1;
    }
    else {
      if (verbose)
	std::cout << cNode->data << "setting value 2 " << std::endl;
      value = 2;
    }
    return branchFound;
  }
}
