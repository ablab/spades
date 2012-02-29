/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "paths/OffsetTracker.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"
#include <queue>

template class SmallVec< ImplicitOffset, MempoolAllocator<ImplicitOffset> >;
template class OuterVec<ImplicitOffsetVec>;

MutableOffsetTracker::MutableOffsetTracker( const vecUnipathSeq& unipathSeqs, 
                                            const MuxGraph& inverseMuxGraph,
                                            const int firstSuperSeq,
                                            const vec<int>& maxOffsets )
{
  cout << "Finding offsets among superseqs in no more than " << flush;
  int maxOffset = Max( maxOffsets );
  
  this->resize( unipathSeqs.size() );
  priority_queue<ExplicitOffset,vector<ExplicitOffset>,
                 ExplicitOffset::OrderByDecreasingAmount> offsetQueue;
  
  for ( vecUnipathSeq::size_type i = firstSuperSeq; i < unipathSeqs.size(); ++i )
    if ( ! unipathSeqs[i].empty() ) {
      this->Add( i, i, 0 );
      offsetQueue.push( ExplicitOffset( i, i, 0 ) );
    }
  
  int kmersPerDot = 100;
  cout << maxOffset/kmersPerDot + 1 << " passes:" << endl;
  
  int nextDot = 0;
  
  while ( ! offsetQueue.empty() ) {
    // Propogate smallest offset.
    ExplicitOffset smallestOffset = offsetQueue.top();
    offsetQueue.pop();
    
    int from = smallestOffset.GetFrom();
    int oldTo = smallestOffset.GetTo();
    int oldAmount = smallestOffset.GetAmount();
    
    while ( oldAmount >= nextDot ) {
      Dot( cout, nextDot/kmersPerDot );
      nextDot += kmersPerDot;
    }
    
    // Where to push it to?
    vec<Mux> inverseMuxes;
    inverseMuxGraph.GetMuxesOf( OrientedKmerPathId( oldTo, false ), inverseMuxes );
    
    for ( vec<Mux>::iterator iMux = inverseMuxes.begin(); iMux != inverseMuxes.end(); ++iMux ) {
      int newTo = iMux->GetPathId().GetId();
      int newAmount = oldAmount + iMux->GetNumKmers();
      
      if ( newAmount > maxOffsets[from] )
        continue;
      
      if ( this->Add( from, newTo, newAmount ) )
        offsetQueue.push( ExplicitOffset( from, newTo, newAmount ) );
    }
  }
  cout << endl;
  
  cout << "done." << endl;
}  


OffsetTracker::OffsetTracker( const MutableOffsetTracker& mutableTracker )
{
  this->ConvertFrom( mutableTracker );
}

void
OffsetTracker::ConvertFrom( const MutableOffsetTracker& mutableTracker ) 
{
  int totalSize = 0;
  for ( int i = 0; i < mutableTracker.size(); ++i )
    totalSize += mutableTracker.GetOffsetsTo(i).size();
  
  m_offsets.Reserve( totalSize, mutableTracker.size() );
  for ( int i = 0; i < mutableTracker.size(); ++i ) {
    const set<ImplicitOffset>& otherOffsets = mutableTracker.GetOffsetsTo(i);
    ImplicitOffsetVec theseOffsets( otherOffsets.size() );
    copy( otherOffsets.begin(), otherOffsets.end(),
          theseOffsets.begin() );
    m_offsets.push_back( theseOffsets );
  }
}
