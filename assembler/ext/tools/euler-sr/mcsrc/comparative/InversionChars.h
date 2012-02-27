/***************************************************************************
 * Title:          InversionChars.h 
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


#include <set>

#include "CharTree.h"
#include "tree/NewettTree.h"
#include "utils.h"


ssize_t ReadCharsFile(std::string &fileName,
		  std::vector<std::string> &species,
		  ssize_t &numChars,
		  std::vector<ssize_t> &invPositions,
		  std::vector<ssize_t> &invEndPositions,
		  std::vector<ssize_t *> &chars,
		  std::vector<std::vector<std::string> *> &invSpecies );
