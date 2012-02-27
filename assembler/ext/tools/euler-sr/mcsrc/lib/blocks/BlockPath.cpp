/***************************************************************************
 * Title:          BlockPath.cpp 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#include "BlockPath.h"



std::ostream& operator<<(std::ostream &out, const BlockPath &path) {
  ssize_t i;
  std::cout << "path of size " << path.vertices.size() << std::endl;
  for (i = 0; i < path.vertices.size(); i++) {
    out << *path.vertices[i] << std::endl;
  }
  return out;
}
