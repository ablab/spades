/***************************************************************************
 * Title:          EmbossAlign.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef EMBOSS_ALIGN
#define EMBOSS_ALIGN

#include <string>

#include "DNASequence.h"
#include "EmbossAlignment.h"
#include "SRSPairParser.h"

ssize_t EmbossAlign(std::string &commandLine, 
		DNASequence &refSeq, DNASequence &qrySeq, 
		EmbossAlignment &alignment);

#endif
