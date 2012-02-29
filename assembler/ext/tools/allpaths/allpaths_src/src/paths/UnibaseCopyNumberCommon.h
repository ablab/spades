///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file UnibaseCopyNumberCommon.h
 * \author tsharpe
 * \date Dec 16, 2011
 *
 * \brief
 */
#ifndef PATHS_UNIBASECOPYNUMBERCOMMON_H_
#define PATHS_UNIBASECOPYNUMBERCOMMON_H_

#include "Basevector.h"
#include "Vec.h"
#include "system/Types.h"
#include <ostream>

void GetMinLength( const vecbvec& unibases, int& min_length, int& nlongest );
void ComputeBias(
     // INPUTS:
     const int K,
     double occCnPloidy,
     const vec<longlong>& biasOccs,
     const vec<longlong>& biasInst,
     // OUTPUT:
     vec<double>& biasCurveLoc,
     // LOGGING:
     const Bool VERBOSE,
     std::ostream& logout );

#endif /* PATHS_UNIBASECOPYNUMBERCOMMON_H_ */
