///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PERFECT_COUNT_H
#define PERFECT_COUNT_H

#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"

// PerfectCount.  For a given set of sequences, return the number of them that
// are subsumed and perfectly match a given "genome", as defined by a lookup file.
// Sequences having length less than K (as defined in lookup table) are ignored.
// The direction argument determines whether we look only for forward
// alignments, or for alignments in either orientation.
//
// PerfectMark.  Same but return vec<Bool> telling which reads are perfect.

int PerfectCount( const vecbasevector& query, const String& lookup_file, 
     const AlignDir direction );

void PerfectMark( const vecbasevector& query, const String& lookup_file,
     const AlignDir direction, vec<Bool>& perfect );

// PerfectCountPlaces.  For each of a given set of sequences, count the number of
// perfect end-to-end placements on a given "genome".
//
// Inaccurate.  This will double-count placements occurring in overlapping segments
// of the lookup table.

void PerfectCountPlaces( const vecbasevector& query, const String& lookup_file, 
     const AlignDir direction, vec<int>& places );

/// PerfectPick.  For each read, pick at random one of its perfect placements on
/// a "genome".  Return a vector consisting of all the places on the genome that
/// are so picked.
///
/// Note: This will slightly overrepresent placements that occur twice because they
/// lie in overlaps between two segments of the lookup table.
///
/// There is a another version in which the read ids are also returned
/// (as places[i].second).

void PerfectPick( const vecbasevector& query, const String& lookup_file,
		  const AlignDir direction, vec<placement_mark>& places,
		  vec<Bool> & queryHasPerfect );

inline void PerfectPick( const vecbasevector& query, const String& lookup_file,
		  const AlignDir direction, vec<placement_mark>& places )
{    vec<Bool> queryHasPerfect;
     PerfectPick( query, lookup_file, direction, places, queryHasPerfect );    }

void PerfectPick( const vecbasevector& query, const String& lookup_file,
		  const AlignDir direction, vec< pair<placement_mark,int> >& places,
		  vec<Bool> & queryHasPerfect );

#endif
