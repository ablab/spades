/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/UnipathSeqBuilder.h"


UnipathSeqBuilder::UnipathSeqBuilder( const vecKmerPath* pUnipaths,
                                      const KmerPathDatabase* pUnipathDB,
                                      ostream* pLog )
  : m_pUnipaths( pUnipaths ),
    m_pUnipathDB( pUnipathDB ),
    m_pLog( pLog )
{
}


void UnipathSeqBuilder::Build( const KmerPath& path,
                               UnipathSeq& theSeq,
                               Mux& theMux ) const
{
  theSeq.clear();
  theMux = Mux();

  if ( path.IsEmpty() )
    return;
  
  KmerPathLoc pathLoc = path.Begin();
  
  const vecKmerPath& unipaths = *m_pUnipaths;
  const KmerPathDatabase& unipathDB = *m_pUnipathDB;

  vec<longlong> answer;
  while ( 1 )
  {
    // Find the unipath that contains the current kmer in the path.
    unipathDB.Contains( pathLoc.GetKmer(), answer );
    
    // Each kmer in the reads should appear in the unipathDB exactly once.
    ForceAssertEq( answer.size(), 1u );
      
    const tagged_rpint& entry = unipathDB[ answer.front() ];
    int unipathId = entry.PathId();
    
    theSeq.push_back( unipathId );
    
    const KmerPath& unipath = unipaths[ unipathId ];
    
    // If we haven't found the left overhang yet, then this must be
    // the first unipath, in which case we should find out how many
    // more kmers there are to the left of the start of the current
    // path.  This is used to form subsumption records later.
    if ( theMux.GetSegment() == -1 ) {
      KmerPathLoc pathStartOnUnipath( unipath, entry.PathPos() );
      pathStartOnUnipath.SetKmer( pathLoc.GetKmer() );
      
      OrientedKmerPathId okpid( 0, false );
      int segment = entry.PathPos();
      int leftOverhang = pathStartOnUnipath - unipath.Begin();
      
      theMux = Mux( okpid, segment, leftOverhang );
    }

    /** 
     * The following code has been replaced by what follows, which
     * should have the identical functionality while being a bit
     * faster.
     
     KmerPathLoc unipathLoc( unipath, entry.PathPos() );
     unipathLoc.SetKmer( pathLoc.GetKmer() );
     
     bool isMatch = IsPerfectMatchToRight( pathLoc, unipathLoc );
     ForceAssert( isMatch );
     
     * If unipaths have been constructed correctly, we can assume
     * that the path and unipath match to the end, so we can just
     * move the pathLoc that far without actually checking that
     * they match, as IsPerfectMatchToRight() does.
     */
    
    KmerPathLoc unipathEnd = unipath.End();
      
    // How many segments to the end of the unipath from here?
    int numSegsToMove = unipathEnd.GetIndex() - entry.PathPos();
    
    // If the unipath ends at or before the last segment of the
    // path, jump there; otherwise, the unipath extends past the
    // end of the path and we're done.
    if ( numSegsToMove > 0 )
      if ( pathLoc.GetIndex() + numSegsToMove < path.NSegments() )
        pathLoc.SetIndex( pathLoc.GetIndex() + numSegsToMove );
      else
        break;
    
    // If the unipath's last kmer is inside the current segment of
    // the path, jump there; otherwise, we're in the last segment
    // of the path and the unipath extends past the end of the
    // path and we're done.
    if ( unipathEnd.GetKmer() <= pathLoc.Stop() )
      pathLoc.SetKmer( unipathEnd.GetKmer() );
    else
      break;
    
    // If the unipath and the path end at the same kmer, we're done.
    if ( pathLoc == path.End() )
      break;
    
    // Otherwise, jump to the next kmer in the path and start again.
    do {
      pathLoc.Increment();
    } while ( pathLoc.isGap() );
  }
}


void UnipathSeqBuilder::Build( const vecKmerPath& kmerPaths,
                               vecUnipathSeq& unipathSeqs,
                               vec<Mux>& muxes ) const
{
  int numPaths = kmerPaths.size();

  int pathsPerDot = 1;
  while ( numPaths / pathsPerDot > 100 )
    pathsPerDot *= 10;

  int numDots = ( numPaths + pathsPerDot - 1 ) / pathsPerDot;

  if ( m_pLog )
    *m_pLog << "Processing paths in " << numDots << " passes." << endl;

  for ( size_t pathId = 0; pathId < kmerPaths.size(); ++pathId )
  {
    if ( m_pLog && pathId % pathsPerDot == 0 )
      Dot( *m_pLog, pathId / pathsPerDot );
    
    UnipathSeq theSeq;
    Mux theMux;
    this->Build( kmerPaths[pathId], theSeq, theMux );

    unipathSeqs.push_back_reserve( theSeq, 0, 2.0 );
    muxes.push_back( theMux );
  }

  if ( m_pLog )
    *m_pLog << endl;
}  
  
