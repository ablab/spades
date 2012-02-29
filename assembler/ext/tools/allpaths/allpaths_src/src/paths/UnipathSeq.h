/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_UNIPATHSEQ_H
#define PATHS_UNIPATHSEQ_H

#include "paths/KmerPath.h"
#include "Intvector.h"

typedef IntVec UnipathSeq;
typedef VecIntVec vecUnipathSeq;

inline
ostream& operator<< ( ostream& out, const UnipathSeq& unipathSeq ) {
  for ( UnipathSeq::size_type i = 0; i < unipathSeq.size(); ++i ) {
    if ( i > 0 ) 
      out << ".";
    out << BaseAlpha( unipathSeq[i] );
  }
  return out;
}

// Given a set of unipaths, convert the given UnipathSeq to a KmerPath.
void ConvertUnipathSeqToKmerPath( const UnipathSeq& unipathSeq,
                                  const vecKmerPath& unipaths,
                                  KmerPath& kmerPath );

// Given a set of unipaths, convert the given UnipathSeq to a KmerPath.
void ConvertUnipathSeqToKmerPath( const UnipathSeq& unipathSeq,
                                  const vec<KmerPath>& unipaths,
                                  KmerPath& kmerPath );

// Given a set of unipaths, convert the given set of UnipathSeqs to a set of KmerPaths.
void ConvertUnipathSeqsToKmerPaths( const vecUnipathSeq& unipathSeqs,
                                    const vecKmerPath& unipaths,
                                    vecKmerPath& kmerPaths );

// Given a set of unipaths, convert the given set of UnipathSeqs to a set of KmerPaths.
void ConvertUnipathSeqsToKmerPaths( const vecUnipathSeq& unipathSeqs,
                                    const vec<KmerPath>& unipaths,
                                    vecKmerPath& kmerPaths );

#endif
