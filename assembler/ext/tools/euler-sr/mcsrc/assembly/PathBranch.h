/***************************************************************************
 * Title:          PathBranch.h 
 * Author:         Mark Chaisson
 * Created:        2008
 * Last modified:  09/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef PATH_BRANCH_H_
#define PATH_BRANCH_H_

#include <stdlib.h>
#include <map>

class PathBranch {
 public:
	typedef std::map<ssize_t, PathBranch*> BranchMap;
	BranchMap branches;
	ssize_t count;
	PathBranch() {
		count = 0;
	}
};

void DeletePathTree(PathBranch &branch);
#endif
