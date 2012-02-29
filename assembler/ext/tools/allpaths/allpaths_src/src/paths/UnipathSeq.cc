/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/UnipathSeq.h"

template <class KmerPathVector>
void ConvertOne( const UnipathSeq& unipathSeq,
                 const KmerPathVector& unipaths,
                 KmerPath& kmerPath ) 
{
  kmerPath.Clear();
  for ( UnipathSeq::size_type i = 0; i < unipathSeq.size(); ++i ) {
    kmerPath.Append( unipaths[ unipathSeq[i] ] );
  }
}

template <class KmerPathVector>
void ConvertMany( const vecUnipathSeq& unipathSeqs,
                  const KmerPathVector& unipaths,
                  vecKmerPath& kmerPaths )
{
  kmerPaths.clear();
  longlong kmerPathsRawsize = 0;
  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqs.size(); ++i )
    for ( UnipathSeq::size_type j = 0; j < unipathSeqs[i].size(); ++j )
      kmerPathsRawsize += unipaths[ unipathSeqs[i][j] ].NSegments();
  kmerPaths.Reserve( kmerPathsRawsize, unipathSeqs.size() );
  
  KmerPath kmerPath;
  for ( vecUnipathSeq::size_type i = 0; i < unipathSeqs.size(); ++i ) {
    ConvertOne( unipathSeqs[i], unipaths, kmerPath );
    kmerPaths.push_back( kmerPath );
  }
}


void ConvertUnipathSeqToKmerPath( const UnipathSeq& x, const vecKmerPath& y, KmerPath& z )
{
  ConvertOne( x, y, z );
}

void ConvertUnipathSeqsToKmerPaths( const UnipathSeq& x, const vec<KmerPath>& y, KmerPath& z )
{
  ConvertOne( x, y, z );
}

void ConvertUnipathSeqsToKmerPaths( const vecUnipathSeq& x, const vecKmerPath& y, vecKmerPath& z )
{
  ConvertMany( x, y, z );
}

void ConvertUnipathSeqsToKmerPaths( const vecUnipathSeq& x, const vec<KmerPath>& y, vecKmerPath& z)
{
  ConvertMany( x, y, z );
}
