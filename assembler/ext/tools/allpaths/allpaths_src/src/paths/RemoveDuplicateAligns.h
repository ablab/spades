///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REMOVE_DUPLICATE_ALIGNS_H
#define REMOVE_DUPLICATE_ALIGNS_H

#include "PairsManager.h"
#include "SeqInterval.h"
#include "VecTemplate.h"
#include "paths/Alignlet.h"

/**
 * RemoveDuplicateAligns
 *
 * Alignments of duplicate molecules are removed by resetting the
 * index of the alignments to -2 ("no alignment found").
 */
void RemoveDuplicateAligns( const PairsManager &pairs,
			    const vec<alignlet> &aligns,
			    vec<int> &index,
			    ostream &log );

#endif
