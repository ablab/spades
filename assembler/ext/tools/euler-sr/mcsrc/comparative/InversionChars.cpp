/***************************************************************************
 * Title:          InversionChars.cpp 
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


#include "InversionChars.h"


ssize_t ReadCharsFile(std::string &fileName,
		  std::vector<std::string> &species,
		  ssize_t &numChars,
		  std::vector<ssize_t> &invStartPositions,
		  std::vector<ssize_t> &invEndPositions,
		  std::vector<ssize_t *> &chars,
		  std::vector<std::vector<std::string> *> &invSpecies ) {

  std::ifstream in;

  openck(fileName, in);

  ssize_t numSpecies, numInvSpec;
  in >> numSpecies >> numChars;

  ssize_t i, start, end, count, number;
  for (i = 0; i < numChars; i++) {
    in >> number >> start >> end >> count;
    invStartPositions.push_back(start);
    invEndPositions.push_back(end);
  }
  
  std::string speciesName;
  ssize_t j, ch;
  ssize_t *posVect;
  std::vector<std::string> *invSpec;

  for (i = 0; i < numSpecies; i++) {
    in >> speciesName;
    species.push_back(speciesName);
    posVect   = new ssize_t[numChars];
    for (j = 0; j < numChars; j++)
      in >> posVect[j];
    chars.push_back(posVect);
  }

  for (i = 0; i < numChars; i++ ) {
    invSpec = new std::vector<std::string>;
    in >> numInvSpec;
    for (j = 0; j < numInvSpec; j++) {
      in >> speciesName;
      invSpec->push_back(speciesName);
    }
    invSpecies.push_back(invSpec);
  }
}
