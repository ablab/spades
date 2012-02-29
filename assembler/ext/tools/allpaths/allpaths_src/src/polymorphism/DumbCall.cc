/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "polymorphism/DumbCall.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"

void SumDumbCalls(VecDumbcallVec& A, VecDumbcallVec& B)
{
    if ( A.size() < B.size() )
        A.resize(B.size());

    for ( unsigned long i = 0; i < A.size(); i++ )
    {
        if ( A[i].size() < B[i].size() )
            A[i].resize(B[i].size());

        for ( unsigned int j = 0; j < A[i].size(); j++ )
            A[i][j] += B[i][j];
    }
}

template class SmallVec< dumbcall, MempoolAllocator<dumbcall> >;
template class OuterVec<DumbcallVec>;
