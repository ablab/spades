/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef FIRST_LOOKUP_H
#define FIRST_LOOKUP_H

#include "lookup/LookAlign.h"
#include "lookup/LookupTable.h"
#include "lookup/FirstLookupFilter.h"

/**
 * FirstLookup
 *
 * Align (possibly chimeric) reads to a reference. Code stolen from
 * BluntJumpAlign.
 *
 * The user should create a FirstLookupFilter, which speeds up FirstLookup
 * by filtering alignments.
 * If no FirstLookupFilter struct is explicitly supplied, the default value is
 * an empty FirstLookupFilter, which does no filtering.
 *
 * aligns: output
 */
void FirstLookup( const vecbasevector& query,
		  const String& lookup_file,
		  vec<look_align>& aligns,
		  const FirstLookupFilter & filter = FirstLookupFilter( ),
		  const int NUM_THREADS = 1);


#endif
