/***************************************************************************
 * Title:          BlastResult.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef BLAST_RESULT_H_
#define BLAST_RESULT_H_


#include "BlastHSP.h"
#include <vector>

class BlastQueryMatch {
public:
  std::string queryName;
  std::vector<BlastHSP> hsps;
};

typedef std::vector<BlastQueryMatch*> BlastResult;

#endif
