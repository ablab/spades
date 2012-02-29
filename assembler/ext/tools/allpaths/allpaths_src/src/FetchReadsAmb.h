/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef FETCHREADSAMB
#define FETCHREADSAMB

#include "Bitvector.h"
#include "CoreTools.h"

enum AMB_STYLE {AMB_EQ_ANY_AMBIGUOUS, AMB_EQ_n, AMB_EQ_Nn, AMB_EQ_LOWER_CASE};

void FetchReadsAmb( vecbitvector& b, String fasta_file, 
     AMB_STYLE amb_style = AMB_EQ_ANY_AMBIGUOUS );

#endif
