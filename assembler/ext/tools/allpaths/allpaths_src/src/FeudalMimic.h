///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Mimic: given a vecbasevector, generate a vecbitbector which has the same
// "dimensions", so that there is a bijective correspondence between bases in one
// and bits in the other.  (Also set the bitvector entries to False.)  Other similar
// feudal conversions go here.
// See FeudalMimic.cc for some instantiations of these templatized functions
// with specific object types.

#ifndef FEUDAL_MIMIC_H
#define FEUDAL_MIMIC_H

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"

// Generic Mimic function.
template<class vecvec1, class vecvec2>
void Mimic( const vecvec1 & in, vecvec2 & out ) {
  out.resize( in.size( ) );
  for (typename vecvec1::size_type i=0; i != in.size( ); ++i)
    out[i].resize(in[i].size( ));
}

// Mimic the shape, but reserve space on the inner vectors, rather than resizing
template<class vecvec1, class vecvec2>
void MimicReserve( const vecvec1 & in, vecvec2 & out )
{
    typedef typename vecvec1::const_iterator I1;
    typedef typename vecvec2::iterator I2;
    out.resize(in.size());
    I2 i2(out.begin());
    for ( I1 i1(in.begin()), iE(in.end()); i1 != iE; ++i1, ++i2 )
    {
        i2->clear();
        i2->reserve(i1->size());
    }
}

// Generic Mimic function, with templatized initialization value.
template<class vecvec1, class vecvec2, class number>
void Mimic( const vecvec1 & in, vecvec2 & out, number value ) {
  out.resize( in.size( ) );
  for (typename vecvec1::size_type i=0; i != in.size( ); ++i)
    out[i].resize(in[i].size( ),value);
}

#endif
