///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SCAFFOLD_ACCURACY_CORE
#define SCAFFOLD_ACCURACY_CORE

#include "Basevector.h"
#include "Bitvector.h"
#include "CoreTools.h"

void SelectTestSequences( int& SAMPLE, const int D, const vecbasevector& assembly, 
     const vecbitvector& assembly_amb, vecbasevector& query,
     vec<int>& t1s, vec<int>& t2s, vec<int>& start1s, vec<int>& start2s );

#endif
