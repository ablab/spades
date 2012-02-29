
///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef LONG_READ_PATCHER_OPTIMIZER_H
#define LONG_READ_PATCHER_OPTIMIZER_H


#include "paths/LongReadTools.h"


void patcher_optimal(const BaseVecVec & unibases,
                     const vec<GapPatcher> & patchers,
                     const size_t ip_best, 
                     const int L,
                     const int nb_rad_SW,
                     const int sz_padding_min,
                     const unsigned i_gap,
                     GapPatcher0 * patcher0_opt_p,
                     const unsigned PATCH_VERBOSITY,
                     vec<double> * timers_p);

#endif
