/***************************************************************************
 * Title:          NetReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  09/06/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef NET_READER_H_
#define NET_READER_H_

#include "Net.h"
#include "NetFile.h"

#include <iostream>

class NetReader {
  static ssize_t ParseLevel(std::ifstream &in, ssize_t &level);
  static void ParseKeywordValues(std::ifstream &in, Net &net);
public:
  static void ReadNetFile(std::string &inFileName, NetFile &netFile);
};

#endif
