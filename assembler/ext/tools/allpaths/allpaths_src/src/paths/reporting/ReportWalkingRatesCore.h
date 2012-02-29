///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__REPORTING__REPORT_WALKING_RATES_CORE_H
#define PATHS__REPORTING__REPORT_WALKING_RATES_CORE_H

#include "Basevector.h"
#include "TokenizeString.h"
#include "VecUtilities.h"

/**
 * ReportWalkingRatesCore
 * 
 * VERBOSE: verbose log flag
 * ids: vector of seed ids
 * info: number of pairs of inserts (walked, total) for each nhood
 */
void ReportWalkingRatesCore( ostream &out,
			     const bool VERBOSE,
			     const vec<int> &ids,
			     const vec< pair<int,int> > &info);

#endif
