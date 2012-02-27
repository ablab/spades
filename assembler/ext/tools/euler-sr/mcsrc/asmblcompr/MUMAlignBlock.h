/***************************************************************************
 * Title:          MUMAlignBlock.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef MUM_ALIGN_BLOCK_H_
#define MUM_ALIGN_BLOCK_H_

#include <fstream>
#include <iostream>

class MUMAlignBlock {
public:
  ssize_t refStrand, qryStrand;
  std::vector<ssize_t> refPos, qryPos, length, refGap, qryGap;
  ssize_t size() { return refPos.size();}
};


#endif
