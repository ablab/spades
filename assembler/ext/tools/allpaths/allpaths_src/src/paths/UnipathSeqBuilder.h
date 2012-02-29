/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_UNIPATHSEQBUILDER_H
#define PATHS_UNIPATHSEQBUILDER_H

#include "Vec.h"

#include "paths/KmerPath.h"
#include "paths/KmerPathDatabase.h"
#include "paths/Mux.h"
#include "paths/UnipathSeq.h"

class UnipathSeqBuilder {
 public:
  UnipathSeqBuilder( const vecKmerPath* pUnipaths, 
                     const KmerPathDatabase* pUnipathDB,
                     ostream* pLog = 0 );

  void Build( const KmerPath& kmerPath,
              UnipathSeq& unipathSeq,
              Mux& theMux ) const;

  void Build( const vecKmerPath& kmerPaths,
              vecUnipathSeq& unipathSeqs,
              vec<Mux>& muxes ) const;

 private:
  const vecKmerPath* m_pUnipaths;
  const KmerPathDatabase* m_pUnipathDB;
  ostream* m_pLog;
};

#endif
