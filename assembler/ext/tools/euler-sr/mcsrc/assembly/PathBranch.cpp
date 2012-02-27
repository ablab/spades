/***************************************************************************
 * Title:          PathBranch.cpp 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "PathBranch.h"


void DeletePathTree(PathBranch &branch) {
	// Perform a post-order traversal of the tree 
	// to delete the leaves.
	PathBranch::BranchMap::iterator branchIt;
	for (branchIt = branch.branches.begin(); 
			 branchIt != branch.branches.end();
			 ++branchIt) {
		if ((*branchIt).second != NULL) {
			DeletePathTree(*(*branchIt).second);
			delete (*branchIt).second;
		}
	}
	branch.branches.clear();
}


