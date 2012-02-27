/***************************************************************************
 * Title:          InversionDB.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  01/08/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef INVERSION_DB_H
#define INVERSION_DB_H

#include <set>
#include <vector>
#include <string>

#include "blocks/BlockDB.h"
#include "ValidatedInversion.h"

void VerifyInversion(DBQuery &query,
		     std::string refSpecies,
		     std::string qrySpecies,
		     std::string sequence,
		     StringSet &temporaryTables, ssize_t maxTempTables,
		     ssize_t startPos, ssize_t endPos,
		     ssize_t doLocal,
		     ssize_t &type,
		     std::map<std::string, DNASequence*> &seqMap,
		     FloatMatrix &scoreMat);

void BinInversionList(DBQuery &query, 
		      InversionList &invList, 
		      std::string seqName,
		      std::vector<ssize_t> &startPos,
		      std::vector<ssize_t> &endPos,
		      BinMap &binnedInversions);

#endif
