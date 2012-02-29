///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Charvector.cc
 * \author tsharpe
 * \date Sep 3, 2009
 *
 * \brief
 */
#include "Charvector.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

template class SmallVec< char, MempoolAllocator<char> >;
template class OuterVec<charvector>;
template class SmallVec< unsigned char, MempoolAllocator<unsigned char> >;
template class OuterVec<ucharvector>;

void StripNewlines( const charvector &in, charvector &out )
{
    out.reserve(in.size()).clear();
    charvector::const_iterator end(in.end());
    for ( charvector::const_iterator itr(in.begin()); itr != end; ++itr )
    {
        char val = *itr;
        if ( val != '\n' )
            out.push_back(val);
    }
}
