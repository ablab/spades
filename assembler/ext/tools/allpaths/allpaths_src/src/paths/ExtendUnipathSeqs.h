/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_EXTENDUNIPATHSEQS
#define PATHS_EXTENDUNIPATHSEQS

#include "paths/KmerPath.h"
#include "paths/Mux.h"
#include "paths/UnipathSeq.h"

/// Given a set of unipaths and UnipathSeqs, extend each element of the
/// given set of UnipathSeqs where the extension is unambiguous.  If,
/// for example, every UnipathSeq that contains some unipath A either
/// has unipath A as the last unipath or has unipath B as the next
/// unipath, then all of the unipaths that end in A may safely be
/// extended to include B.
///
/// Note: It is safe to pass the same vecUnipathSeq object for both the
/// unipathSeqs and extendedSeqs arguments.
void ExtendUnipathSeqs( const vecKmerPath& unipaths,
                        const vecUnipathSeq& unipathSeqs, 
                        vecUnipathSeq& extendedSeqs,
                        vec<Mux>& muxes );

/// This version allows you to specify a different set of unipath seqs
/// to extend from the total collection of unipath seqs.
void ExtendUnipathSeqs( const vecKmerPath& unipaths,
                        const vecUnipathSeq& unipathSeqs,
                        const vecUnipathSeq& unipathSeqsToExtend, 
                        vecUnipathSeq& extendedSeqs,
                        vec<Mux>& muxes );
#endif
