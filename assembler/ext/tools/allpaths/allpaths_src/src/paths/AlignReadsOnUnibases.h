///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__ALIGN_READS_ON_UNIBASES_H
#define PATHS__ALIGN_READS_ON_UNIBASES_H

#include "Basevector.h"

void AlignReadsOnUnibases( // input
			   const vecbvec &jbases,
			   const vecbvec &unibases,
			   
			   // output
			   vec< triple<int64_t,int,int> > &jaligns,
			   vec<basevector> &jbases_sorted,
			   vec<int64_t> &jbases_sorted_id,
			   
			   // log stream
			   ostream &out );

#endif
