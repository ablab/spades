/***************************************************************************
 * Title:          BlockPath.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _BLOCK_PATH_
#define _BLOCK_PATH_

#include <iostream>
#include "Block.h"

// Use BlockPath for storing a path in the graph.  I'm not sure yet if
// this will be useful.

class BlockPath {
public:
  std::vector<Block*> vertices;  // The blocks themselves
  std::vector<ssize_t>    indices; // The index of this sequence in the block
  ssize_t size() { assert(vertices.size() == indices.size()); return vertices.size(); }
}; 

std::ostream& operator<<(std::ostream &out, const BlockPath &path);

typedef std::vector<Block*>::iterator PathIterator;
typedef std::vector<ssize_t>::iterator IndexIterator;
#endif
