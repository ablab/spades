///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file Basevector.cc
 * \author tsharpe
 * \date Sep 23, 2009
 *
 * \brief
 */
#include "Basevector.h"

void ReverseComplement( vecbasevector& vbv )
{
    vecbvec::iterator end(vbv.end());
    for ( vecbvec::iterator itr(vbv.begin()); itr != end; ++itr )
        itr->ReverseComplement();
}

#include "feudal/OuterVecDefs.h"
template class OuterVec<bvec>;
template class OuterVec<bvec,bvec::allocator_type>;
template class OuterVec< OuterVec<bvec,bvec::allocator_type>,
                         MempoolOwner<unsigned char> >;

