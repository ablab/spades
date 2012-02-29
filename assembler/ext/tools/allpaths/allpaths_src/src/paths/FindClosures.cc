/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "ReadPairing.h"
#include "paths/AddSuperReads.h"
#include "paths/FindClosures.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/KmerPathMuxSearcher.h"
#include "paths/PairedPair.h"

void FindClosures( const vec<pp_pair>& ppp,
                   const vec<Bool>& pairs_to_close,
                   const vec< vec<pp_closure> >& prior_closures,
                   const HyperKmerPath& h,
                   const double sdMult,
                   vec< HyperKmerPath >& closures,
                   vec<Bool>& fail,
                   const unsigned int max_pseudo_closures,
                   const unsigned int max_closures )
{
  fail.resize_and_set( ppp.size( ), False );
  vec<int> L( h.EdgeObjectCount( ) );
  for ( int i = 0; i < L.isize( ); i++ )
    L[i] = h.EdgeLength(i);
  vecKmerPath paths;
  vec<read_pairing> pairs;
  for ( int i = 0; i < ppp.isize( ); i++ ) {
    const pp_pair& p = ppp[i];
    read_pairing rp;
    rp.t = other;
    rp.weight = 1;
    rp.id1 = paths.size( );
    static KmerPath x;
    x.Clear( );
    for ( int j = 0; j < p.LeftSize( ); j++ )
      x.Append( h.EdgeObject( p.Left(j) ) );
    paths.push_back_reserve(x, 0, 2.0);
    double oldgap = p.Gap( ), olddev = p.Dev( ), m = sdMult;
    int newgap = int( round(oldgap) ), newdev = int( ceil(olddev) );
    while(1) {
      Bool raise = False;
      if ( !( oldgap - m * olddev >= newgap - int( floor( m * newdev ) ) ) )
        raise = True;
      if ( !( oldgap + m * olddev <= newgap + int( ceil( m * newdev ) ) ) )
        raise = True;
      if (raise)
        ++newdev;
      else
        break;
    }
    if ( newdev == 0 ) ++newdev;
    rp.sep = newgap - h.K( ) + 1;
    rp.sd = newdev;
    rp.id2 = paths.size( );
    x.Clear( );
    for ( int j = 0; j < p.RightSize( ); j++ )
      x.Append( h.EdgeObject( p.Right(j) ) );
    paths.push_back_reserve(x, 0, 2.0);
    if ( pairs_to_close[i] ) 
      pairs.push_back(rp);    
  }
  cout << "Reserving" << endl;
  int pathsSize = paths.size();
  for ( size_t i = 0; i < prior_closures.size(); ++i )
    pathsSize += prior_closures[i].size();

  paths.reserve( pathsSize );
  cout << "Pushing back" << endl;
  for ( int i = 0; i < prior_closures.isize(); ++i )
    for ( int j = 0; j < prior_closures[i].isize(); ++j ) {
      const pp_closure& p = prior_closures[i][j];
      static KmerPath x;
      x.Clear();
      for ( int k = 0; k < p.isize(); ++k )
        x.Append( h.EdgeObject( p[k] ) );
      paths.push_back(x);
    }
  cout << "Calling FindClosures()" << endl;
  FindClosures( paths, pairs, sdMult, h.K( ), closures, fail, 
                max_pseudo_closures, max_closures );
}


void PrepareData( const vecKmerPath& paths,
                  const vec<read_pairing>& pairs,
                  const double sdMult,
                  const int K,
                  vec<read_pairing>& newPairs,
                  vecKmerPath& allPathsFw,
                  vecKmerPath& allPathsRc,
                  MuxGraph& allMuxes,
                  SubsumptionList& allSubs,
                  OffsetTracker* pOffTracker )
{
  const int numInserts = pairs.size();
  
  vecKmerPath pathsFw;
  vecKmerPath pathsRc;
  
  vec<bool> isPaired( paths.size(), false );

  newPairs.reserve( numInserts );

  longlong pathsFwRawsize = 0;
  longlong pathsRcRawsize = 0;
  int pathsSize = numInserts;
  for ( int pass = 0; pass < 2; ++pass ) {
    if ( pass == 1 ) {
      pathsFw.Reserve( pathsFwRawsize, pathsSize );
      pathsRc.Reserve( pathsRcRawsize, pathsSize );
    }
    for ( int i = 0; i < numInserts; ++i ) {
      const read_pairing& origPair = pairs[i];
      if ( pass == 0 ) {
        pathsFwRawsize += paths[ origPair.id1 ].NSegments();
        pathsRcRawsize += paths[ origPair.id2 ].NSegments();
        isPaired[ origPair.id1 ] = true;
        isPaired[ origPair.id2 ] = true;
      }
      else {
        read_pairing newPair = origPair;
        newPair.id1 = pathsFw.size(), newPair.id2 = pathsRc.size();
        pathsFw.push_back( paths[ origPair.id1 ] );
        pathsRc.push_back( paths[ origPair.id2 ] );
        newPairs.push_back( newPair );
      }
    }
    if ( pass == 0 ) {
      for ( size_t i = 0; i < paths.size(); ++i )
        if ( ! isPaired[i] ) {
          pathsFwRawsize += paths[i].NSegments();
          ++pathsSize;
        }
    } else {
      for ( size_t i = 0; i < paths.size(); ++i )
        if ( ! isPaired[i] ) {
          pathsFw.push_back( paths[i] );
          pathsRc.push_back( KmerPath() );
        }
    }
  }

  const int maxSep = 0; // no max sep

  cout << "Calling AddSuperReads" << endl;
  AddSuperReads( pathsFw, pathsRc, newPairs, K, maxSep, Float(sdMult),
                 allPathsFw, allPathsRc, allMuxes, allSubs, pOffTracker );
}


void FindClosures( const vecKmerPath& paths,
                   const vec<read_pairing>& pairs,
                   const double sdMult,
                   const int K,
                   vec< HyperKmerPath >& closures, 
                   vec<Bool>& fail,
                   const unsigned int max_pseudo_closures,
                   const unsigned int max_closures )

{
  fail.resize_and_set( pairs.size( ), False );
  const int numInserts = pairs.size();
  
  vecKmerPath pathsFw;
  vecKmerPath pathsRc;
  vec<read_pairing> newPairs;
  vecKmerPath allPathsFw, allPathsRc;
  MuxGraph allMuxes;
  SubsumptionList allSubs;

  cout << "Calling PrepareData" << endl;
  PrepareData( paths, pairs, sdMult, K,
               newPairs, allPathsFw, allPathsRc, allMuxes, allSubs, 0 );

  ReadFillDatabase identityFillDb;
  vec<int> readLengths( allPathsFw.size() );
  for ( size_t i = 0; i < allPathsFw.size(); ++i )
    readLengths[i] = allPathsFw[i].KmerCount();
  
  int verbosity = 0;
  KmerPathMuxSearcher searcher( &allMuxes, &identityFillDb, &allSubs, &readLengths, verbosity );

  cout << "Walking " << numInserts << " inserts:" << endl;
  for ( int i = 0; i < numInserts; ++i ) {
    cout << i << ":";
    const read_pairing& thePair = newPairs[i];
    
    int expLength = readLengths[ thePair.id1 ] + thePair.sep + (K-1);
    int minLength = expLength - (int)ceil( sdMult * (double)thePair.sd );
    int maxLength = expLength + (int)ceil( sdMult * (double)thePair.sd );

    cout << " walking..." << flush;
    MuxSearchResult result;
    searcher.FindClosures( thePair.id2, thePair.id1,
                           minLength, maxLength,
                           result );

    cout << " calculating..." << flush;
    
    if (result.hit_search_limit) 
      fail[i] = True;   
    if ( max_pseudo_closures > 0 && result.num_closures_found > (int) max_pseudo_closures )
      fail[i] = True;

    closures.push_back( HyperKmerPath( K, vec<KmerPath>() ) );
    if ( ! fail[i] ) {
      result.WalkGraph( ).MakeHyperKmerPath( &allPathsFw, &allPathsRc, &allSubs, closures.back() );
    }

    cout << " done." << endl;
  }
}


void FindClosureLengths( const vecKmerPath& paths,
                         const vec<read_pairing>& pairs,
                         const double sdMult,
                         const int K,
                         vec< vec<int> >& closureLengths )
{
  const int numInserts = pairs.size();
  
  vecKmerPath pathsFw;
  vecKmerPath pathsRc;
  vec<read_pairing> newPairs;
  vecKmerPath allPathsFw, allPathsRc;
  MuxGraph allMuxes;
  SubsumptionList allSubs;
  OffsetTracker offTracker;

  PrepareData( paths, pairs, sdMult, K,
               newPairs, allPathsFw, allPathsRc, allMuxes, allSubs, &offTracker );

  closureLengths.clear();
  closureLengths.resize( numInserts );

  for ( int i = 0; i < numInserts; ++i ) {
    int id1 = newPairs[i].id1;
    int id2 = newPairs[i].id2;
    
    vec<Mux> muxes;
    allMuxes.GetMuxesOf( OrientedKmerPathId( id1, false ), muxes );
    ForceAssertEq( muxes.size(), 1u );
    Mux leftMux = muxes.front();
    int leftSuperId = leftMux.GetPathId().GetId();
    int leftInset = leftMux.GetNumKmers();

    allMuxes.GetMuxesOf( OrientedKmerPathId( id2, true ), muxes );
    ForceAssertEq( muxes.size(), 1u );
    Mux rightMux = muxes.front();
    int rightSuperId = rightMux.GetPathId().GetId();
    int rightOffset = rightMux.GetNumKmers() + allPathsRc[id2].KmerCount() + (K-1);

    OffsetTracker::Range range;
    range = offTracker.GetOffsetsToFrom( rightSuperId, leftSuperId );

    for ( ; range.first != range.second; ++range.first ) {
      closureLengths[i].push_back( range.first->GetAmount() - leftInset + rightOffset );
    }
  }
}
