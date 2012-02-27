/***************************************************************************
 * Title:          AxtLib.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef _AXTLIB_
#define _AXTLIB_
#include <assert.h>
#include "AxtEntry.h"
#include "lav/LAVFile.h"

void AxtToLAV(AxtEntries &axtList, LAVFile &lavFile, ssize_t refSeqLen, ssize_t qrySeqLen);  
void PrintAxt(AxtEntries &axtFile, std::ostream &out);
#endif
