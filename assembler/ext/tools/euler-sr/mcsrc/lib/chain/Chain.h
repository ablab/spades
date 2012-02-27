/***************************************************************************
 * Title:          Chain.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/26/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CHAIN_H_
#define CHAIN_H_

#include "ChainHeader.h"
#include <vector>
class Chain {
public:
  ChainHeader  header;
  std::vector<ssize_t> size, dt, dq;
  ssize_t id;
  ssize_t numAlign() {
    return size.size();
  }
  Chain & operator=(const Chain &src) {
		if (this != &src) {
			header.score = src.header.score;
			header.tName = src.header.tName;
			header.tSize = src.header.tSize;
			header.tStrand = src.header.tStrand;
			header.tStart  = src.header.tStart;
			header.tEnd    = src.header.tEnd;
			header.qName   = src.header.qName;
			header.qSize   = src.header.qSize;
			header.qStrand = src.header.qStrand;
			header.qStart  = src.header.qStart;
			header.qEnd    = src.header.qEnd;
			id = src.id;
			// TODO: why aren't other fields copied too? size, dt, dq
		}
    return *this;
  }
};

#endif
