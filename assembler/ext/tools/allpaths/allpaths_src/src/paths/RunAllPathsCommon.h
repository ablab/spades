///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Contains additional code for RunAllPathsLG.
*/

#ifndef RUN_ALLPATHS_COMMON_H
#define RUN_ALLPATHS_COMMON_H

#include "String.h"


// Create or update genome.size file in reference directory
void UpdateRefGenomeSize(String ref_dir);


// Sets threading values for each module, returns global thread count
int ParseThreadArgs(const String& THREADS,
		    int& LR_THREADS, int& CP_THREADS, int& FF_THREADS,
		    int& ECJ_THREADS, int& KP_THREADS,
		    int& CUG_THREADS, int& RDR_THREADS, int& FE_THREADS,
		    int& MN_THREADS, int& PR_THREADS, int& PC_THREADS);


void DisplaySystemStatus();
  
#endif
