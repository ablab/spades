/***************************************************************************
 * Title:          ChainReader.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CHAIN_READER_H_
#define CHAIN_READER_H_

#include <string>
#include <vector>
#include <iostream>
#include <fstream>

#include "Chain.h"
#include "ChainHeader.h"


class ChainReader {
  static ssize_t ReadChainHeader(std::ifstream &in, Chain &chain);
  static ssize_t ReadAlignment(std::ifstream &in, Chain &chain);
public:
  static void ReadChainFile(std::string chainFileName, 
			   std::vector<Chain*> &chains);
};


#endif
