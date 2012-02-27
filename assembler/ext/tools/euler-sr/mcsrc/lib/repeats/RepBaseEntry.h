/***************************************************************************
 * Title:          RepBaseEntry.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef REP_BASE_ENTRY_H_
#define REP_BASE_ENTRY_H_

#include <string>
#include "DNASequence.h"

class RepBaseEntry {
public:
  std::string name;
  std::string repClass;
  DNASequence reference;
  ssize_t ParseRepBaseEntry(std::ifstream &in);
};

#endif
