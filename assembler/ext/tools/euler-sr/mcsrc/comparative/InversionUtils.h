/***************************************************************************
 * Title:          InversionUtils.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  11/27/2009
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef INVERSION_UTILS_H_
#define INVERSION_UTILS_H_
#include <iostream>
#include <fstream>
#include <string>
#include <iterator>
#include <set>

// 3rd party
#include "mysql/mysql.h"

// Mine 
#include "net/NetDB.h"
#include "blocks/BlockDB.h"
#include "ValidatedInversion.h"
#include "Inversion.h"
#include "InversionChars.h"
#include "InversionFile.h"
#include "emboss/EmbossAlign.h"
#include "DNASequence.h"
#include "SeqReader.h"
#include "TupleLib.h"
#include "blocks/dbOrthoPos.h"
#include "lav/LAVTable.h"

#define WITH_INVERSION 0
#define WITHOUT_INVERSION 1
#define UNKNOWN 2
#define DELETED 3

void DoLocalAlign(DNASequence &refSeq, DNASequence &qrySeq,
		  ssize_t invStart, ssize_t  invEnd, ssize_t qryStart, ssize_t qryEnd,
		  ssize_t &result);

ssize_t CheckForDeletion(ssize_t tStart, ssize_t tEnd, ssize_t qStart, ssize_t qEnd, double maxRat);


void ReadValidatedFile(std::string &fileName, StringVector &speceis, InversionList &invList);


#endif
