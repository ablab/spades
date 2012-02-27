/***************************************************************************
 * Title:          ChainUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef CHAIN_UTILS_H_
#define CHAIN_UTILS_H_

#include "chain/Chain.h"
#include "lav/LAVFile.h"

void ChainToAlignment(Chain &chain, LAVAlignedContig &alignment);
void ChainsToAlignFile(std::vector<Chain*> &chains, LAVFile &lavFile);



#endif
