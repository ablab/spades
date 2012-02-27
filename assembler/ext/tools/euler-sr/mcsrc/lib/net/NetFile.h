/***************************************************************************
 * Title:          NetFile.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NETFILE_H_
#define NETFILE_H_

#include <vector>
#include <string>

#include "Net.h"

class NetFile {
public:
  std::string name;
  ssize_t length;
  std::vector<Net*> nets;
  ssize_t size() { return nets.size(); }
};


#endif
