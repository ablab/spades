///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SUCK_SCAFFOLDS_H
#define SUCK_SCAFFOLDS_H

#include "PairsManager.h"
#include "Superb.h"
#include "paths/Alignlet.h"
#include "paths/HyperFastavector.h"

/**
 * SuckScaffolds
 *
 * Absorb small scaffolds based only on linking evidence (provided the
 * small scaffolds fit inside an existing gap of some other scaffold).
 * Returns the number of sucking events.
 *
 * pairs: pairs info
 * aligns_index: index from read it to align object (or <0 if not unique)
 * aligns: alignments of reads onto contigs (these may be flipped)
 * scaffolds: both input and output
 * rctig: stores info on which contigs are rc-ed (scaffolds can be rc-ed)
 * plog: optional log stream
 */
int SuckScaffolds( const PairsManager &pairs,
		   const vec<int> &aligns_index,
		   vec<alignlet> &aligns,
		   vec<superb> &scaffolds,
		   vec<Bool> &rctig,
		   ostream *plog = 0 );

#endif
