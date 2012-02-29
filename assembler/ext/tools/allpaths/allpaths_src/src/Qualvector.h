///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/// \file
/// This file defines the typedef "qualvector", which stores quality scores
/// between 0 and 255 as a vector of unsigned chars, and the typedef
/// "vecqualvector", which stores a vector of qualvectors.
/// \ingroup grp_quals

#ifndef QUALVECTOR
#define QUALVECTOR

#include "feudal/SerfVec.h"
#include "feudal/MasterVec.h"
#include "String.h"
#include <ostream>

/// Logical type for quality scores
typedef uint8_t qual_t;

/// Vector of quality scores, for example representing the quality of each base
/// in one read.
typedef SerfVec<qual_t> QualVec;
typedef QualVec qualvector;
typedef QualVec qvec;


/// Vector of vectors of quality scores, for example representing the quality
/// of each base in each read in a set of reads.
typedef MasterVec< SerfVec<qual_t> > QualVecVec;
typedef QualVecVec vecqualvector;
typedef QualVecVec vecqvec;

typedef OuterVec< OuterVec<qvec,MempoolAllocator<qual_t> >,
                  MempoolOwner<qual_t> > qvec3;

///Produces fasta format quals, mirrors basevector::Print()
void Print( std::ostream &out, const qualvector &q, const String &name,
            const int scores_per_line = 25 );

/// CopyQuals: copy elements from one qualvector to another. If rev_from=True,
/// then copy from reverse(from), starting at from_start on reverse(from).
/// This mirrors CopyBases in Basevector.h.

inline void CopyQuals( const qualvector& from, int from_start, qualvector& to,
     int to_start, int count, bool rev_from = false )
{    if ( !rev_from )
     {    for ( int i = 0; i < count; i++ )
               to[ to_start + i ] = from[ from_start + i ];    }
     else
     {    for ( int i = 0; i < count; i++ )
               to[ to_start + i ]
                    = from[ from.size( ) - (from_start + i) - 1 ];    }    }

// ReadFastaQuals: read quality scores from a fasta quality score file.  If
// ids_to_read is supplied, it should be a sorted list of indices of the records
// to be read.

template <class T> class vec;

void ReadFastaQuals( const String& fn, vecqualvector& qual,
                        const vec<int>* ids_to_read = 0 );

#endif
