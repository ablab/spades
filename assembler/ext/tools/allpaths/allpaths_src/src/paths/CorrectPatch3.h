///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef CORRECT_PATCH3_H
#define CORRECT_PATCH3_H

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"
#include "paths/AssemblyEdit.h"

Bool CorrectPatch3( const basevector& LEFT, const basevector& RIGHT,
     const vecbasevector& fbases, const vecqualvector& fquals, 
     const PairsManager& fpairs, const vec< kmer<20> >& fheads, 
     const vec<int64_t>& fids, assembly_edit& e, ostringstream& out,
     const Bool verbose, const vecbasevector& genome2, const int LG,
     const vec< vec< pair<int,int> > > & Glocs );

#endif
