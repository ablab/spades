///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SMITHWATBANDEDA
#define SMITHWATBANDEDA

#include "Alignment.h"
#include "Basevector.h"

// Warning.  An object of type X contains the error count at a given position in
// the dynamic programming matrix.  Below, we use X = unsigned char, which could
// easily overflow.  Use of a larger type is recommended but carries with it a
// performance penalty.

// Note also the extremely ugly inclusion of the dummy argument x, which is needed,
// for unknown reasons.

template<class X> float SmithWatBandedA2( const basevector& S, 
     const basevector& T, int offset, int bandwidth, align& a, int& errors, 
     ostream* log = 0, int mis=2, int ins=3, int del=3 );

inline float SmithWatBandedA( const basevector& S, const basevector& T, 
     int offset, int bandwidth, align& a, int& errors, ostream* log = 0, 
     int mis=2, int gap=3  )
{    
     return SmithWatBandedA2<unsigned char>( 
          S, T, offset, bandwidth, a, errors, log, mis, gap, gap  );    }

inline float SmithWatBandedA( const basevector& S, const basevector& T, 
     int offset, int bandwidth, alignment& a, ostream *log = 0 )
{    align temp; int errors;
     float result = SmithWatBandedA( S, T, offset, bandwidth, temp, errors, log );
     a.Set( temp, errors );
     return result;    }

#endif
