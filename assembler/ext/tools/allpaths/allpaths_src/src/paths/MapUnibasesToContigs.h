///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__MAP_UNIBASES_TO_CONTIGS__H
#define PATHS__MAP_UNIBASES_TO_CONTIGS__H

#include "Basevector.h"
#include "Intvector.h"
#include "util/SearchFastb2Core.h"

/**
 * MapUnibasesToContigs
 *
 * Generate a unibases_to_contigs map, using SearchFastb2Core to align
 * contigs to unibases. Temp files are deleted at the end.
 *
 * OUTBASE: output is saved as <OUTBASE>.{u2c} (as a UInt64VecVec)
 * log: used to log progress
 */
void MapUnibasesToContigs( const int K,
			   const String &unibases_file,
			   const String &contigs_file,
			   const String &OUTBASE,
			   ostream *log = 0 );

#endif
