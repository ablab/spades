///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef FETCHREADS
#define FETCHREADS

#include "Basevector.h"
#include "Floatvector.h"
#include "CoreTools.h"
#include "Qualvector.h"

/// \file FetchReads.h
/// Routines for getting reads from a fasta file into mastervecs.
/// (vecbasevector and vecString).
/// Can also use FastFetchReads in FastaFileset.h, which is twice as fast.

/// Put reads from a fasta file into a vecbasevector and a names vector.
/// Ambiguous bases are saved as random bases.
/// If the names vector pointer is empty, ignores names and thus preserves
/// efficiency.
///
/// The vector b is cleared first.

void FetchReads( vecbasevector& b, vecqualvector& q, vecString * names,
                 unsigned int n, String fasta_file,
                 int amb_break = 0, int min_size = 0,
                 ostream& out = cout,
                 Bool no_q = False, const vec<int>* ids_to_read = 0 );

inline void FetchReads( vecbasevector& b, vecString& names, String fasta_file )
{    vecqualvector q;
     FetchReads( b, q, &names, 0, fasta_file, 0, 0, cout, True );    }

///Put reads from a fasta file into a vecbasevector.
/// Ambiguous bases are saved as random bases.
void FetchReads( vecbasevector& b, vecqualvector& q, unsigned int n,
                 String fasta_file, int amb_break = 0, int min_size = 0,
                 ostream& out = cout,
                 Bool no_q = False, const vec<int>* ids_to_read = 0 );

///Put reads from a fasta file into a vecbasevector.
/// Ambiguous bases are saved as random bases.
void FetchReads( vecbasevector& b, unsigned int n, String fasta_file,
     int amb_break = 0, int min_size = 0, ostream& out = cout,
     const vec<int>* ids_to_read = 0 );

/// Read a "fasta" file whose entries are floats.

void FetchReads( VecFloatVec& f, const String& fasta_file );

#endif
