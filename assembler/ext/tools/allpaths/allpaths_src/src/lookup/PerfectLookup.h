///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PERFECT_LOOKUP_H
#define PERFECT_LOOKUP_H

#include "Basevector.h"
#include "CoreTools.h"
#include "lookup/LookAlign.h"

/**
   \file
   
// Given a vector "query" of sequences, find all perfect alignments
// to target (defined by "lookup_file", from MakeLookupTable with the given "K"), 
// where the default behaviour is that the alignment is required to go end-to-end 
// on the query.  Queries of length < K are not aligned. The option "direction"
// specifies whether only forward alignments are to be found, or if reverse
// alignments are to be found too.  Output is "aligns".
//
// If the flag subsumed_only is set to False, then partial alignments will also be
// found using the same low kmer freq method described below. This does NOT return
// all possible partial alignments, just those 'seeded' by the low frequency kmer
// positions in the target. However, if a set of overlapping target basevectors
// are used, with a fixed overlap (in bases) given in target_seq_overlap, then all
// the partial alignments for each overlap region WILL be found. For example, if 
// the target was created from a hyperbasevector, then all possible alignments of
// the query sequences to the hyperbasevector will be found. The value of the
// target_seq_overlap in this case would be K-1 (wherre K is for the associated
// hyperkmerpath, not the lookup table).
//
// Note that for certain applications, the code be speeded up by keeping the
// lookup table in memory, rather than reading it from disk on each call.
//
// Method: for each query sequence, we find a kmer in it that appears a minimal
// number of times in the target.  Each occurrence of this kmer in the target 
// defines a possible offset, which we test directly.
//
// I'm not sure how this handles Ns in the target.
//
// WARNING: this code will return duplicate alignments.

   \sa ImperfectLookup()
   \ingroup grp_eval
*/

/**
   Direction of alignment of a query to a target.
   \sa PerfectLookup()
*/
enum AlignDir {
  FW, /**< find only forward alignments */
  FW_OR_RC  /**< find both forward alignment and reverse complement alignments */
};

/// \copydoc PerfectLookup.h
/// Aligns all queries to the target

void PerfectLookup( const unsigned int K, const vecbasevector& query, 
		    const String& lookup_file, vec<look_align>& aligns, 
                    const AlignDir direction, const Bool subsumed_only = True, 
                    const unsigned int target_seq_overlap = 0 );

/// \copydoc PerfectLookup.h
/// Aligns only queries in the range specified by range_start and range_end
/// Note that [range_start,range_end] is a CLOSED interval.

void PerfectLookup( const unsigned int K, const vecbasevector& query, 
		    const String& lookup_file, vec<look_align>& aligns,
		    const AlignDir direction, int range_start, int range_end,
		    const Bool subsumed_only = True, 
                    const unsigned int target_seq_overlap = 0 );

// Parallelized versions of PerfectLookup.  Something is wrong with them because
// they will occasionally fail.

void ParallelPerfectLookup( const int max_threads, const unsigned int K, 
     const vecbasevector& query, const String& lookup_file, vec<look_align>& aligns, 
     const AlignDir direction, const Bool subsumed_only = True, 
     const unsigned int target_seq_overlap = 0 );

void ParallelPerfectLookup( const int max_threads, const unsigned int K, 
     const vecbasevector& query, const String& lookup_file, vec<look_align>& aligns,
     const AlignDir direction, int range_start, int range_end,
     const Bool subsumed_only = True, const unsigned int target_seq_overlap = 0 );

#endif
