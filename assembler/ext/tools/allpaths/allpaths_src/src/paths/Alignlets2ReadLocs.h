///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef ALIGNLETS_2_READ_LOCS_H
#define ALIGNLETS_2_READ_LOCS_H

#include "PairsManager.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "paths/Alignlet.h"
#include "paths/ReadLoc.h"

/**
 * Alignlets2ReadLocs
 *
 * Convert pairs, alignlets, and matching index file into a (sorted!)
 * vector of read_locs. WARNING! All locs will be classified as type 1
 * (jump).
 */
void Alignlets2ReadLocs( const PairsManager &pairs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 vec<read_loc> &locs,
			 ostream *log = 0 );

/**
 * LoadReadLocs
 *
 * Load vectors of alignlets and indexes, and convert them into a
 * vector of read_locs. See WARNING in Alignlets2ReadLocs.
 *
 * pairs_file: the PairsManager file
 * head_file: it needs <head>.{,index}
 * locs: output
 * log: optional log stream (to show progress)
 */
void LoadReadLocs( const String &pairs_file,
		   const String &head_file,
		   vec<read_loc> &locs,
		   ostream *log = 0 );

#endif
