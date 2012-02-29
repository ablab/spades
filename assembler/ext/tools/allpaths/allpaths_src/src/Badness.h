///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


#include <math.h>

#include "system/Assert.h"

// ================================================================================
//
// Badness(overlap, errors): return a quality measure for two reads with given
// putative overlap and error count
//
// ================================================================================

inline float Badness( int overlap, int errors )
{    Assert( overlap > 0 );
     return 1000.0 * float(1+errors)/pow(float(overlap), float(1.5));    }
