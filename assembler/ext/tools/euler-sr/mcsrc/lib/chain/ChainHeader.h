/***************************************************************************
 * Title:          ChainHeader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CHAIN_HEADER
#define CHAIN_HEADER

#include <string>

class ChainHeader {
public:
  ssize_t score;
  std::string tName;
  ssize_t tSize, tStrand;
  ssize_t tStart, tEnd;
  std::string qName;
  ssize_t qSize, qStrand;
  ssize_t qStart, qEnd;
};


#endif
