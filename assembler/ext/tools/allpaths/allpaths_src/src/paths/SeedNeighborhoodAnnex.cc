///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


/*******************************************************************************
 *
 *        SEED NEIGHBORHOOD ANNEX
 *
 * This module contains the definitions for several member functions in the
 * SeedNeighborhood class.  Much of the work done in the SeedNeighborhood is
 * performed by calls to these functions.  Most of these functions are versions
 * of already-existing functions that have been adapted for the RunAllPathsLG
 * pipeline.
 *
 * The file SeedNeighborhoodAnnex.h contains no functions.
 * For function declarations, see SeedNeighborhood.h.
 *
 *
 * Josh Burton
 * September 2009
 *
 ******************************************************************************/

#include <set>
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "paths/GetNexts.h"
#include "paths/HyperBasevector.h"
#include "paths/KmerPath.h"
#include "paths/InternalMerge.h"
#include "paths/InsertWalker.h"
#include "paths/PdfEntry.h"
#include "paths/SeedNeighborhood.h"
#include "paths/SimpleLoop.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h" // BuildUnipathAdjacencyHyperKmerPath
#include "paths/UnipathNhoodCommon.h" // fsepdev, ustart











Bool Linked( longlong x, longlong y, const vec<ReadLocationLG>& ulocs,
     const vec<longlong>& uindex, const PairsManager& pairs )
{    vec<longlong> yreads;
     yreads.clear( );
     longlong indx1 = uindex[x], indx2 = uindex[x+1];
     longlong indy1 = uindex[y], indy2 = uindex[y+1];
     for ( longlong i = indy1; i < indy2; i++ )
          yreads.push_back( ulocs[i].ReadId( ) );
     for ( longlong j1 = indx1; j1 < indx2; j1++ )
     {    const ReadLocationLG& rl1 = ulocs[j1];
          if ( !rl1.Fw( ) ) continue;
          longlong id1 = rl1.ReadId( );
	  longlong pi = pairs.getPairID( id1 );
	  if ( pi < 0 ) continue;
	  longlong id2 = pairs.getPartner( pi, id1 );
          int j2 = BinPosition( yreads, id2 );
          if ( j2 < 0 ) continue;
          const ReadLocationLG& rl2 = ulocs[ indy1 + j2 ];
          if ( !rl2.Rc( ) ) continue;
          return True;    }
     return False;    }








/*******************************************************************************
 *
 * OBJECT MEMBER FUNCTIONS (declared in SeedNeighborhood.h)
 *
 ******************************************************************************/





// Given a seed unipath, find all nearby <normal unipaths> to form a
// <neighborhood> around this seed.
void
SeedNeighborhood::FindUnipathCloud( const digraphE<fsepdev> & LG,
				    const vec<int> & predicted_CNs,
				    const vec<Bool> & branch,
				    const int MAX_DEV )
{
  Bool verbose = _NHOOD_VERBOSITY & VERBOSITY_UNIPATH_CLOUD;
  const int MIN_KMERS = 25; // ignore unipaths shorter than this
  const int MIN_BRANCH_KMERS = 100; // ignore branch unipaths shorter than this
  const int MAX_TO_CUT = 4000; // cutpoints longer than this are kept
  const int MAX_PATHS_IN_NHOOD = 400; // how big the nhood can get
  
  // Require the seed to be sufficiently long.
  ForceAssertGe( (*_global_unipath_lengths)[_seed_ID], MIN_KMERS );
  
  vec<int> zero_dev( 1, 0 );
  vec<ustart> unprocessed;
  vec<int> allowed;
  set<int> announced;
  
  if (verbose) {
    announced.insert( _seed_ID );
    cout << "starting with seed unipath " << _seed_ID << endl;
  }

  for ( int upass = 1; upass <= 2; upass++ ) {
    _cloud_unipaths.clear( ), unprocessed.clear( );
    unprocessed.push_back( ustart( _seed_ID, 0, zero_dev ) );
    while( unprocessed.nonempty( ) ) {
      pop_heap( unprocessed.begin(), unprocessed.end(), ustart::OrderByDescendingMeanDev() );
      int w = unprocessed.back( ).Uid( );
      int wstart = unprocessed.back( ).Start( );
      vec<int> wdev = unprocessed.back( ).Dev( );
      // If we've processed this unipath before, that ustart had a
      // MeanDev no worse than this one, so skip this one.
      Bool used = False;
      for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
        if ( w == _cloud_unipaths[j].Uid( ) ) {
          used = True;
          break;
        }
      if ( used ) {
        unprocessed.pop_back();
        continue;
      }

      _cloud_unipaths.push_back( unprocessed.back( ) );
      unprocessed.pop_back();

      if ( predicted_CNs[w] > Max(_mcn_other, predicted_CNs[_seed_ID]) )
           continue;

      for ( int i = 0; i < LG.From(w).isize( ); i++ ) {
        int x = LG.From(w)[i];
	if ( (*_global_unipath_lengths)[x] < MIN_KMERS ) continue;
        if ( predicted_CNs[_seed_ID] == _ploidy && _mcn_other <= _ploidy 
             && branch[x] && (*_global_unipath_lengths)[x] <= MIN_BRANCH_KMERS )
        {    continue;    }
        if ( upass == 2 && !BinMember( allowed, x ) ) continue;
        int wxsep = LG.EdgeObjectByIndexFrom( w, i ).Sep( );
        int xstart = wstart + (*_global_unipath_lengths)[w] + wxsep;
        if ( xstart > (*_global_unipath_lengths)[_seed_ID] + _NHOOD_RADIUS ) continue;

        vec<int> dev = wdev;
        dev.push_back( LG.EdgeObjectByIndexFrom( w, i ).Dev( ) );
        ustart new_ustart( x, xstart, dev );
        if ( new_ustart.MeanDev( ) > MAX_DEV ) continue;

        Bool used = False;
        for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
          if ( x == _cloud_unipaths[j].Uid( ) ) {
            used = True;
            break;
          }
        if (used) continue;

        if ( verbose && !Member( announced, x ) ) 
        {    announced.insert(x);
             cout << "accepting unipath " << x << " via link from " << w << endl;   }

        unprocessed.push_back( new_ustart );
        push_heap( unprocessed.begin(), unprocessed.end(), 
             ustart::OrderByDescendingMeanDev() );
        
        if ( _cloud_unipaths.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;
      }
      if ( _cloud_unipaths.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;

        for ( int i = 0; i < LG.To(w).isize( ); i++ ) {
          int x = LG.To(w)[i];
	  if ( (*_global_unipath_lengths)[x] < MIN_KMERS ) continue;
          if ( predicted_CNs[_seed_ID] == 1 && _mcn_other <= 1 
               && branch[x] && (*_global_unipath_lengths)[x] <= MIN_BRANCH_KMERS )
          {    continue;    }
          if ( upass == 2 && !BinMember( allowed, x ) ) continue;
          int xwsep = LG.EdgeObjectByIndexTo( w, i ).Sep( );
          int xstart = wstart + - xwsep - (*_global_unipath_lengths)[x];
          if ( xstart + (*_global_unipath_lengths)[x] < -_NHOOD_RADIUS ) continue;

          vec<int> dev = wdev;
          dev.push_back( LG.EdgeObjectByIndexTo( w, i ).Dev( ) );
          ustart new_ustart( x, xstart, dev );
          if ( new_ustart.MeanDev( ) > MAX_DEV ) continue;

          Bool used = False;
          for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
            if ( x == _cloud_unipaths[j].Uid( ) ) {
              used = True;
              break;
            }
          if (used) continue;

          if ( verbose && !Member( announced, x ) ) 
          {    announced.insert(x);
               cout << "accepting unipath " << x << " via link to " << w << endl;   }

          unprocessed.push_back( new_ustart );
          push_heap( unprocessed.begin(), unprocessed.end(), 
               ustart::OrderByDescendingMeanDev() );

          if ( _cloud_unipaths.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) 
               break;
        }
      
      if ( _cloud_unipaths.isize( ) + unprocessed.isize( ) >= MAX_PATHS_IN_NHOOD ) break;
    }

    // If there are any unprocessed unipaths left, that are not already in 
    // processed, add them in.

    for ( int j = 0; j < unprocessed.isize( ); j++ )
    {    Bool used = False;
         for ( int z = 0; z < _cloud_unipaths.isize( ); z++ )
              if ( _cloud_unipaths[z].Uid( ) == unprocessed[j].Uid( ) ) used = True;
         if ( !used ) 
         {    if (verbose) 
                   cout << "adding in unipath " << unprocessed[j].Uid( ) << endl;
              _cloud_unipaths.push_back( unprocessed[j] );    }    }
    Sort(_cloud_unipaths);

    // Note: it is likely that by more effectively using the available 
    // information, we could come up with better estimates for the positions
    // of the unipaths in the neighborhood.
    
    // Screen the neighborhood to remove unipaths which are incompatible
    // with the seed neighborhood (v).  This is an interim step.  We will
    // want to systematically address incompatibilities.  Also remove cut
    // points.
    
    if ( upass == 1 ) {
      vec<ustart> processedx;
      
      for ( int j = 0; j < _cloud_unipaths.isize( ); j++ ) {
        int start = _cloud_unipaths[j].Start( ), w = _cloud_unipaths[j].Uid( );
        int dev = _cloud_unipaths[j].MeanDev( );
        if ( start >= 0 ) {
          Bool linked = Linked( _seed_ID, w, (*_global_unilocs), (*_global_unilocs_index), *_global_pairs );
	  if ( linked )
            processedx.push_back( _cloud_unipaths[j] );    
        }
        else {
          Bool linked = Linked( w, _seed_ID, (*_global_unilocs), (*_global_unilocs_index), *_global_pairs );
	  if ( linked )
            processedx.push_back( _cloud_unipaths[j] );
        }
      }
      _cloud_unipaths = processedx;
      for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
        allowed.push_back( _cloud_unipaths[j].Uid( ) );
      UniqueSort(allowed);
      continue;
    }

    // For each graph cut point corresponding to a unipath of length
    // less than a middling insert length, remove the vertices that do
    // not lie in the connected component of the neighborhood seed v.
    
    if ( upass == 2 ) {
      // Create a graph that contains all the selected unipaths and
      // the edges between them.
      vec<int> processedv;
      
      for ( int j = 0; j < _cloud_unipaths.isize( ); j++ ) {
	if ( (*_global_unipath_lengths)[ int(_cloud_unipaths[j].Uid( )) ] < MIN_KMERS ) continue;
	processedv.push_back( int(_cloud_unipaths[j].Uid( )) );
      }
      digraphE<fsepdev> N( LG, processedv );

      // Find the vertices (other than the seed) that, if removed,
      // would disconnect the graph.
      vec<int> cuts;
      N.CutPoints(cuts);
      vec<Bool> to_remove;
      to_remove.resize_and_set( N.N( ), False );
      int seed_vx = Position( processedv, _seed_ID );
      for ( int j = 0; j < cuts.isize( ); j++ ) {
        int c = cuts[j];
        if ( c == seed_vx ) continue;
        int cut_unipath = processedv[c];
        if ( (*_global_unipath_lengths)[ cut_unipath ] > MAX_TO_CUT ) continue;
        equiv_rel e( N.N( ) );
        for ( int vx = 0; vx < N.N( ); vx++ ) {
          if ( vx == c ) continue;
          for ( int j = 0; j < N.From(vx).isize( ); j++ ) {
            int w = N.From(vx)[j];
            if ( w == c ) continue;
            e.Join(vx, w);
          }
        }
        vec<int> keep;
        e.Orbit( seed_vx, keep );
        Sort(keep);
        for ( int i = 0; i < N.N( ); i++ ) {
          if ( i != c && !BinMember( keep, i ) )
            to_remove[i] = True;
        }
      }
      EraseIf( _cloud_unipaths, to_remove );    

    }
  }
}





// Add in unipaths of predicted copy number <= ploidy
void
SeedNeighborhood::ExpandUnipathCloud( const digraphE<fsepdev> & LG,
				      const vec<int> & predicted_CNs,
				      const int MAX_DEV )
{
     vec<ustart> processed_extra;
     processed_extra.clear( );
     for ( int i = 0; i < _cloud_unipaths.isize( ); i++ )
     {    int w = _cloud_unipaths[i].Uid( );
          if ( predicted_CNs[w] > _ploidy ) continue;
          int wstart = _cloud_unipaths[i].Start( );
          vec<int> wdev = _cloud_unipaths[i].Dev( );
          for ( int i = 0; i < LG.From(w).isize( ); i++ )
          {    int x = LG.From(w)[i];
               if ( predicted_CNs[x] > 2*_ploidy ) continue;
               int wxsep = LG.EdgeObjectByIndexFrom( w, i ).Sep( );
               int xstart = wstart + (*_global_unipath_lengths)[w] + wxsep;
               if ( xstart > (*_global_unipath_lengths)[_seed_ID] + _NHOOD_RADIUS ) continue;
               Bool used = False;
               for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
               {    if ( x == _cloud_unipaths[j].Uid( ) )
                    {    used = True;
                         break;    }    }
               if (used) continue;
               for ( int j = 0; j < processed_extra.isize( ); j++ )
               {    if ( x == processed_extra[j].Uid( ) )
                    {    used = True;
                         break;    }    }
               if (used) continue;
               vec<int> dev = wdev;
               dev.push_back( LG.EdgeObjectByIndexFrom( w, i ).Dev( ) );
               processed_extra.push_back( ustart( x, xstart, dev ) );
               if ( processed_extra.back( ).MeanDev( ) > MAX_DEV )
                    processed_extra.resize( processed_extra.isize( ) - 1 );    }
          for ( int i = 0; i < LG.To(w).isize( ); i++ )
          {    int x = LG.To(w)[i];
               if ( predicted_CNs[x] > 2*_ploidy ) continue;
               int xwsep = LG.EdgeObjectByIndexTo( w, i ).Sep( );
               int xstart = wstart + - xwsep - (*_global_unipath_lengths)[x];
               if ( xstart + (*_global_unipath_lengths)[x] < -_NHOOD_RADIUS ) continue;
               Bool used = False;
               for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
               {    if ( x == _cloud_unipaths[j].Uid( ) )
                    {    used = True;
                         break;    }    }
               if (used) continue;
               for ( int j = 0; j < processed_extra.isize( ); j++ )
               {    if ( x == processed_extra[j].Uid( ) )
                    {    used = True;
                         break;    }    }
               if (used) continue;
               vec<int> dev = wdev;
               dev.push_back( LG.EdgeObjectByIndexTo( w, i ).Dev( ) );
               processed_extra.push_back( ustart( x, xstart, dev ) );
               if ( processed_extra.back( ).MeanDev( ) > MAX_DEV )
                    processed_extra.resize( processed_extra.isize( ) - 1 );    }    }
     
     if ( _NHOOD_VERBOSITY ) {
       cout << "ExpandUnipathCloud is adding the following unipaths:";
       for ( size_t i = 0; i < processed_extra.size(); i++ )
	 cout << "\t" << processed_extra[i].Uid();
       cout << endl;
     }
     for ( int j = 0; j < processed_extra.isize( ); j++ )
          _cloud_unipaths.push_back( processed_extra[j] );
     Sort(_cloud_unipaths);    }







/**
   Function: FindPrimaryReadCloud

   Find the reads which go in a particular neighborhood.   Also return their orientations and predicted positions.

   This function is derived from PopulateNhoodWithReads, in UnipathNhood.cc.

   Rough algorithm:

   We want to identify reads where at least one read of the pair is within <inner radius> of the <neighborhood seed>.
   For each unipath within the <outer radius> of the neighborhood seed, we take reads aligned to that unipath.
   From the estimated distance of the unipath to the nhood seed, and the position of the read on the unipath, we
   estimate how far the read is from the nhood seed; if within the inner radius, accept both the read and its partner.
*/
void
SeedNeighborhood::FindPrimaryReadCloud( const vecbvec& global_reads_bases,
     const Bool& LOCAL_PRIMARY )
{
  const int MAX_DEV = 2000; // how stretchy sep from seed can get
  const double max_read_dist_devs = 2.0;
  _read_IDs.clear( ), _read_locs.clear( ), _pair_IDs.clear( );
  for ( int j = 0; j < _cloud_unipaths.isize( ); j++ ) {
    // take one unipath from the nhood, take its start relative to the seed,
    // and its unipath id.
    int start = _cloud_unipaths[j].Start( );
    int w = _cloud_unipaths[j].Uid( );
    vec<int> dev = _cloud_unipaths[j].Dev( );
    // Look at all the reads aligned to this unipath-from-the-nhood.
    longlong indw1 = (*_global_unilocs_index)[w];
    longlong indw2 = (*_global_unilocs_index)[w+1];
    for ( longlong t = indw1; t < indw2; t++ ) {
      longlong id1 = (*_global_unilocs)[t].ReadId( );
      Bool rc1 = (*_global_unilocs)[t].Rc( );
      // Does this read fall roughly on this unipath?
      int start1 = start + (*_global_unilocs)[t].Start( );
      int stop1 = start1 + (*_global_paths)[id1].KmerCount( ) - 1;
      double dist1 = Distance( start1, stop1, -_NHOOD_RADIUS_INTERNAL, _NHOOD_RADIUS_INTERNAL + (*_global_unipath_lengths)[_seed_ID] );
      double dev1 = _cloud_unipaths[j].MeanDev( );
      Bool used_id1 = False;
      if ( dev1 == 0 || dist1/dev1 <= max_read_dist_devs ) {
	_read_IDs.push_back( make_pair( id1, rc1 ? ORIENT_RC : ORIENT_FW ) );
	_read_locs.push_back(start1);
	used_id1 = True;
      }
      longlong pi = _global_pairs->getPairID( id1 );
      if ( pi < 0 ) continue;
      longlong id2 = _global_pairs->getPartner( pi, id1 );
      int sep12 = _global_pairs->sep(pi), dev12 = _global_pairs->sd(pi);
      double dev2 = double(dev12) * double(dev12);
      for ( int t = 0; t < dev.isize( ); t++ )
	dev2 += double(dev[t]) * double(dev[t]);
      dev2 = sqrt(dev2);
      if ( dev2 > MAX_DEV ) continue;
      int start2;
      if ( !rc1 ) start2 = stop1 + sep12 + _K - 1;
      else start2 = start1 - (sep12 + _K - 1) - (*_global_paths)[id2].KmerCount( );
      int stop2 = start2 + (*_global_paths)[id2].KmerCount( ) - 1;
      double dist2 = Distance( start2, stop2, -_NHOOD_RADIUS_INTERNAL, _NHOOD_RADIUS_INTERNAL + (*_global_unipath_lengths)[_seed_ID] );
      if ( dev2 == 0 || dist2/dev2 <= max_read_dist_devs ) {
	_read_IDs.push_back( make_pair( id2, !rc1 ? ORIENT_RC : ORIENT_FW ) );    
	_read_locs.push_back(start2);
	if (used_id1) _pair_IDs.push_back(pi);
      }
    }
  }
  SortSync(_read_IDs, _read_locs), UniqueSort(_pair_IDs);
  vec<Bool> to_delete( _read_IDs.isize( ), False );
  for ( int i = 0; i < _read_IDs.isize( ) - 1; i++ )
    if ( _read_IDs[i] == _read_IDs[i+1] ) to_delete[i] = True;
  EraseIf( _read_IDs, to_delete ), EraseIf( _read_locs, to_delete );

  if (LOCAL_PRIMARY) {
    ofstream* pout = new ofstream;
    OpenOfstream( *pout, _outhead + ".primary.fasta" );
    for ( int i = 0; i < _read_IDs.isize( ) ; i++ )
      global_reads_bases[_read_IDs[i].first].Print( *pout, "id=" + ToString(_read_IDs[i].first) );
    pout->close();
  }
}







// Find the secondary read cloud ("cloud2"): reads that share a kmer with the
// reads in the primary read cloud ("cloud1").
// This function is adapted from GetShortInsertReads in UnipathNhood.cc.
void
SeedNeighborhood::FindSecondaryReadCloud( )
{
  //cout << "Size of _read_IDs:\t" << _read_IDs.isize( ) << endl;
  
  // 1. Form a list L0 of the KmerPathIntervals which appear in S
  // (where S = _read_IDs).
  
  vec<KmerPathInterval> L0, L;
  for ( int i = 0; i < _read_IDs.isize( ); i++ ) {
    longlong id = _read_IDs[i].first;
    if ( !_read_IDs[i].second )
      for ( int j = 0; j < (*_global_paths)[id].NSegments( ); j++ )
	L0.push_back( (*_global_paths)[id].Segment(j) );
    else
      for ( int j = 0; j < (*_global_paths_rc)[id].NSegments( ); j++ )
	L0.push_back( (*_global_paths_rc)[id].Segment(j) );
  }

  sort( L0.begin( ), L0.end( ), cmp_start );
  //cout << "Size of L0:\t" << L0.isize( ) << endl;
  
  // 2. Condense the list L0 into L.
  
  for ( int i = 0; i < L0.isize( ); i++ ) {
    longlong Istart = L0[i].Start( ), Istop = L0[i].Stop( );
    
    int j;
    for ( j = i + 1; j < L0.isize( ); j++ ) {
      if ( L0[j].Start( ) > Istop ) break;
      Istop = Max( Istop, L0[j].Stop( ) );
    }
    L.push_back( KmerPathInterval( Istart, Istop ) );
    i = j - 1;
  }
  //cout << "Size of L:\t" << L.isize( ) << endl;
  
  // 3. Find all reads R that share kmers with L.
  
  vec< pair<longlong,Bool> > R;
  for ( int i = 0; i < L.isize( ); i++ ) {
    vec<longlong> con;
    Contains( (*_global_pathsdb), L[i], con );
    for ( int j = 0; j < con.isize( ); j++ ) {
      const tagged_rpint& t = (*_global_pathsdb)[ con[j] ];
      longlong id = t.ReadId( );
      Bool rc = ( t.PathId( ) < 0 );
      R.push_back( make_pair( id, rc ) );
    }
  }
  UniqueSort(R);
  //cout << "Size of R:\t" << R.isize( ) << endl;
  
  // 4. For each read in R, work from left to right, trying to cover it
  // by reads in S, yielding a new set P.
  
  vec<tagged_rpint> Spathsdb;
  
  for ( int i = 0; i < _read_IDs.isize( ); i++ ) {
    longlong read_ID = _read_IDs[i].first;
    if ( !_read_IDs[i].second )
      (*_global_paths)   [ read_ID ].AppendToDatabase(Spathsdb, read_ID);
    else
      (*_global_paths_rc)[ read_ID ].AppendToDatabase(Spathsdb, -read_ID-1);
  }
  Prepare(Spathsdb);
  
  vec< pair< longlong, Bool> > cloud2;
  for ( int i = 0; i < R.isize( ); i++ ) {
    longlong id = R[i].first;
    Bool rc = R[i].second;
    const KmerPath& p = ( !rc ? (*_global_paths)[id] : (*_global_paths_rc)[id] );
    if ( SubContig( p, (*_global_paths), (*_global_paths_rc), Spathsdb ) )
      cloud2.push_back( R[i] );
  }
  //cout << "Size of cloud2:\t" << cloud2.isize( ) << endl;
  
  // Merge the secondary read cloud into the first one.  To keep read_IDs and
  // read_locs in sync we add fake entries to the latter, which won't be used.
  _read_IDs.append( cloud2 );
  for ( int i = 0; i < cloud2.isize( ); i++ )
    _read_locs.push_back(0);
  UniqueSortSync( _read_IDs, _read_locs );
}






// An adaptation of SelectPairsToWalk for the LG case.  The algorithm is
// equivalent to calling SelectPairsToWalk (from LocalizeReadsAnnex2)
// with PAIRS_SAMPLE="cover".
// The number of selected pairs will be no greater than (n_pairs_to_select);
// it may be less, if there are not enough logical pairs to choose from.
vec<longlong>
SeedNeighborhood::SelectPairsToWalk( const int n_pairs_to_select, int & n_logical_pairs ) const
{

  // Find the reads in the nhood that meet a cloud unipath.  We use unilocs
  // to answer this question.  It would be better to directly check for overlap
  // of kmers between the reads and the cloud unipaths.

  vec<longlong> placed_ids;
  for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
  {    int w = _cloud_unipaths[j].Uid( );
       longlong indw1 = (*_global_unilocs_index)[w];
       longlong indw2 = (*_global_unilocs_index)[w+1];
       for ( longlong t = indw1; t < indw2; t++ ) 
       {    longlong id1 = (*_global_unilocs)[t].ReadId( );
            placed_ids.push_back(id1);    }    }
  UniqueSort(placed_ids);
  vec< vec<ReadLocationLG> > placements( placed_ids.size( ) );
  for ( int j = 0; j < _cloud_unipaths.isize( ); j++ )
  {    int w = _cloud_unipaths[j].Uid( );
       longlong indw1 = (*_global_unilocs_index)[w];
       longlong indw2 = (*_global_unilocs_index)[w+1];
       for ( longlong t = indw1; t < indw2; t++ ) 
       {    longlong id1 = (*_global_unilocs)[t].ReadId( );
            longlong p1 = BinPosition( placed_ids, id1 );
            placements[p1].push_back( (*_global_unilocs)[t] );    }    }

  // Create form of _cloud_unipaths that is sorted by unipath id.

  vec<int> cloud_uid, cloud_start;
  for ( int i = 0; i < _cloud_unipaths.isize( ); i++ )
  {    cloud_uid.push_back( _cloud_unipaths[i].Uid( ) );
       cloud_start.push_back( _cloud_unipaths[i].Start( ) );    }
  SortSync( cloud_uid, cloud_start );
  
  // Each pair<ho_interval, longlong> is the interval covered by a read pair, then that read pair's ID
  vec<pair<ho_interval, longlong> > covers;
  covers.reserve( _n_pairs );
  //cout << Date() << ": N PAIRS = " << _n_pairs << endl;
  
  // For each local pair, determine if the pair is logical (therefore walkable)
  // and then determine what interval in the local neighborhood it covers
  for ( longlong ptui = 0; ptui < _n_pairs; ++ptui ) {
    longlong pi = _pair_IDs[ptui];
    longlong global1 = _global_pairs->ID1(pi), global2 = _global_pairs->ID2(pi);

    // Ignore pairs for which one or both ends does not meet the unipath cloud.
    int p1 = BinPosition( placed_ids, global1 );
    int p2 = BinPosition( placed_ids, global2 );
    if ( p1 == -1 || p2 == -1 ) continue;
    //cout << "\tPair " << ptui << " passed step 1" << endl;

    int local1 = BinPosition( _local_to_global, global1 );
    int local2 = BinPosition( _local_to_global, global2 );
    ForceAssertGe( local1, 0 );
    ForceAssertGe( local2, 0 );
    
    // If the reads in this pair point in the same direction, ignore the pair
    // Otherwise, re-order them (if necessary) so read1 is FW and read2 is RC
    bool rc1 = _read_IDs[local1].second;
    bool rc2 = _read_IDs[local2].second;
    if ( rc1 == rc2 ) continue;
    //cout << "\tPair " << ptui << " passed step 2" << endl;
    if ( rc1 ) {
      swap( local1, local2 );
      swap( global1, global2 );
    }
    
    int start1 = _read_locs[local1];
    int start2 = _read_locs[local2];
    int len1 = (*_global_read_lengths)[global1];
    int len2 = (*_global_read_lengths)[global2];
    int end1 = start1 + len1;
    int end2 = start2 + len2;
    
    // Ignore pairs in which the reads point away from each other
    // (This is kind of circular.)
    if ( end1 >= start2 ) continue;
    //cout << "\tPair " << ptui << " passed step 3" << endl;

    // Another test for the same thing.
    Bool valid = False;
    for ( int x1 = 0; x1 < placements[p1].isize( ); x1++ )
    {    const ReadLocationLG& rl1 = placements[p1][x1];
         int pu1 = BinPosition( cloud_uid, rl1.Contig( ) );
         int ustart1 = cloud_start[pu1] + rl1.Start( );
         int ustop1 = ustart1 + len1;
         for ( int x2 = 0; x2 < placements[p2].isize( ); x2++ )
         {    const ReadLocationLG& rl2 = placements[p2][x2];
              int pu2 = BinPosition( cloud_uid, rl2.Contig( ) );
              int ustart2 = cloud_start[pu2] + rl2.Start( );
              int ustop2 = ustart2 + len2;
              if ( ustop1 <= ustart2 && rl1.Fw( ) && rl2.Rc( ) ) valid = True;
              if ( ustop2 <= ustart1 && rl2.Fw( ) && rl1.Rc( ) ) 
                   valid = True;    }    }
    if ( !valid ) continue;
    //cout << "\tPair " << ptui << " passed step 4" << endl;
    
    // Ignore pairs that are more than NHOOD_RADIUS_INTERNAL in size
    const int outer = _NHOOD_RADIUS_INTERNAL;
    const int seed_length = (*_global_unipath_lengths)[_seed_ID];
    if ( end1 < -outer || start1 > seed_length + outer ) continue;
    //cout << "\tPair " << ptui << " passed step 5" << endl;
    if ( end2 < -outer || start2 > seed_length + outer ) continue;
    //cout << "\tPair " << ptui << " passed step 6" << endl;
    
    // Record the interval that consists of the space between these reads
    covers.push_back( make_pair( ho_interval( end1, start2 ), ptui ) );
  }
  
  n_logical_pairs = covers.isize( );
  //cout << Date() << ": N COVERS = " << n_logical_pairs << endl;
  if ( n_logical_pairs == 0 ) return vec<longlong>();
  
  // Sort the pairs by ho_interval order
  Sort( covers );
  
  // Choose a total of (n_pairs_to_select) pairs from the set of logical pairs
  // We can choose evenly spaced pairs instead of randomizing, because the
  // pairs are sorted and therefore give reasonably consistent coverage
  vec<longlong> selected_pairs;
  selected_pairs.reserve( n_pairs_to_select );
  
  for (int i = 0; i < n_logical_pairs; i++ ) {
    // This logic guarantees that the first and last pairs will be selected,
    // and that the intermediate pairs will be evenly spaced
    if ( i * (n_pairs_to_select-1) < (n_logical_pairs-1) * selected_pairs.isize() ) continue;
    selected_pairs.push_back( _pair_IDs[ covers[i].second ] );
  }
  //cout << Date() << ": N SELECTED PAIRS = " << selected_pairs.size() << endl;
  
  UniqueSort( selected_pairs );
  return selected_pairs;
}
















// Build an edited version of the HyperKmerPath from which hanging ends and
// cyclic parts have been deleted.
HyperKmerPath
SeedNeighborhood::MakeAcyclicHKP( const vecKmerPath & new_unipaths, const vecbasevector & new_unibases, const int min_size ) const
{
  // Here we start a series of steps that directly build an edited version of
  // the unipath graph from which hanging ends and cyclic parts have been 
  // deleted.  Below this gets merged into the 'final answer'.
  //
  // First build the unipath graph.
     
  vec< vec<int> > from, to( new_unibases.size( ) );
  GetNexts( _K, new_unibases, from );
  for ( size_t u = 0; u < new_unibases.size( ); u++ )
    {    Sort( from[u] ); 
      for ( int j = 0; j < from[u].isize( ); j++ )
	to[ from[u][j] ].push_back(u);    }
  digraph A( from, to );
  HyperKmerPath h;
  BuildUnipathAdjacencyHyperKmerPath( _K, A, new_unipaths, h );

  // Delete edges that are not forward.

  vec<int> to_delete;
  for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
       if ( !_local_unipaths_fw[i] ) to_delete.push_back(i);
  h.DeleteEdges(to_delete);

  // Remove hanging ends.

  RemoveHangingEnds( h, &KmerPath::KmerCount, 250, 5.0 );
  h.RemoveUnneededVertices( );
  h.RemoveSmallComponents(min_size);
  
  // Remove the loop subgraph, and also remove branches as much as possible.
  // This step is the meat of the MakeAcyclicHKP function.
  h.MakeAcyclic();
  
  h.RemoveSmallComponents(min_size);
  h.RemoveDeadEdgeObjects( );

  // Delete reverse complement duplicated components.

  /*
  h.DeleteReverseComplementComponents( );
  h.RemoveDeadEdgeObjects( );
  */
     
  return h;
}



// walk_to_do: Represents an walk between two unipaths that is to be performed.
// Used in WalkInserts, below.
struct walk_to_do {
  int u1, u2; // the unipaths to walk between
  ho_interval range; // the acceptable range of separations from (the start of) u1 to (the start of) u2
  
  // Comparison operators - for sorting
  friend Bool operator==( const walk_to_do & a, const walk_to_do & b ) {
    return ( a.u1 == b.u1 && a.u2 == b.u2 && a.range == b.range );
  }
  
  friend Bool operator<( const walk_to_do & a, const walk_to_do & b ) {
    if ( a.u1 != b.u1 ) return a.u1 < b.u1;
    if ( a.u2 != b.u2 ) return a.u2 < b.u2;
    return a.range < b.range;
  }
};



// WalkInserts: Take a selected set of local read pairs, and walk as many of the
// pairs as possible by traversing the local unipath graph.  Output a HKP
// representing the union of all walked inserts.
//
// Algorithm:
// 1. For each read pair, find the unipaths the reads lie on.
// 2. Sort and condense the list of pairs of unipaths.
// 3. Walk between the pairs of unipaths (InsertWalker::WalkUnipaths).  This is
//    the most time-consuming step by far.
// 4. Combine all of the found walks into a HyperKmerPath.
// 5. Some small steps to clean up the HKP.
HyperKmerPath
SeedNeighborhood::WalkInserts( const vec<longlong> & selected_pair_IDs,
			       const int STRETCH,
			       TaskTimer & timer,
			       int & n_to_walk,
			       int & n_walked,
                               vec<Bool>& walked ) const
{
  HyperKmerPath null_HKP;
  
  // Map of unipaths to their RCs.
  vec<int> to_rc;
  UnipathInvolution( _local_unipaths, _local_unipathsdb, to_rc );
  
  // Set up an InsertWalker!
  InsertWalker insert_walker( _K, &_local_unipaths, &_local_unipathsdb, &_AG );
  
  // List of all unipaths traversed by insert walks, and implied adjacencies.
  // Note that not all adjacencies between unipaths in unipaths_walked will be
  // in unipath_adjs, because some unipaths in unipaths_walked will come from
  // different insert walks.
  vec<Bool> unipaths_walked( _n_local_unipaths, False );
  digraph unipath_adjs( _n_local_unipaths );
  n_walked = 0;
  
  
  // Make a list of walks to do, catalogued by unipath IDs and distance.
  // The walk_to_do struct is defined above.
  vec<walk_to_do> walks_to_do;
  
  // Loop over each of the long-insert pairs that we've selected for walking.
  vec<int> ids;
  walked.resize( selected_pair_IDs.size( ), False );
  for ( int i = 0; i < selected_pair_IDs.isize(); i++ ) {
    longlong pID = selected_pair_IDs[i];
    int sep = _global_pairs->sep( pID ), sd = _global_pairs->sd( pID );
    
    // Find the reads in this pair, expressed as KmerPaths.
    longlong ID1 = BinPosition( _local_to_global, _global_pairs->ID1( pID ) );
    longlong ID2 = BinPosition( _local_to_global, _global_pairs->ID2( pID ) );
    const KmerPath& read1 = _paths   [ID1];
    const KmerPath& read2 = _paths_rc[ID2];
    
    // Find the unipaths landed on by these reads.
    walk_to_do w;
    int sep_offset;
    insert_walker.PlaceReadsOnUnipaths( read1, read2, w.u1, w.u2, sep_offset );
    if ( w.u1 == -1 ) continue;
    // PRINT2( BaseAlpha(w.u1), BaseAlpha(w.u2) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    
    // Find the allowable range of distances from (the start of) u1 to
    // (the start of) u2.
    sep -= sep_offset;
    w.range = ho_interval( sep - sd*STRETCH, sep + sd*STRETCH );
    walks_to_do.push_back( w );
    ids.push_back(i);
  }
  
  
  // Sort and condense the list of walks to do.  This step prevents us from
  // redundantly walking through the same regions over and over again.
  SortSync( walks_to_do, ids );
  vec<Bool> walks_to_erase( walks_to_do.size(), false );
  int last_i = -1, last_u1 = -1, last_u2 = -1;
  
  for ( size_t i = 0; i < walks_to_do.size(); i++ ) {
    const walk_to_do & w = walks_to_do[i];
    
    // Whenever multiple walks_to_do share the same u1 and u2, merge their
    // ranges (even if they don't overlap - which is rare.)
    if ( w.u1 == last_u1 && w.u2 == last_u2 ) {
      walk_to_do & last_w = walks_to_do[last_i];
      last_w.range.SetStart( Min( w.range.Start(), last_w.range.Start() ) );
      last_w.range.SetStop ( Max( w.range.Stop (), last_w.range.Stop () ) );
      walks_to_erase[i] = true;
    }
    else {
      last_i = i;
      last_u1 = w.u1;
      last_u2 = w.u2;
    }
  }
  EraseIf( walks_to_do, walks_to_erase );
  EraseIf( ids, walks_to_erase );
  n_to_walk = walks_to_do.size();
  
  // Keep track of redundant inserts.
  vec<bool> redundant( selected_pair_IDs.size( ), true ); 
  
  // Now, do the walks!
  for ( size_t i = 0; i < walks_to_do.size(); i++ ) {
    const walk_to_do & w = walks_to_do[i];
    redundant[ ids[i] ] = false;
    
    // Walk between these unipaths using the InsertWalker.
    // If the InsertWalker succeeds, it returns the set of unipath IDs walked.
    // This is a runtime bottleneck!
    set<int> this_walk = insert_walker.WalkUnipaths( w.u1, w.u2, w.range, timer, false );

    //if (timer.TimedOut()) return null_HKP;
    
    // Remove any unipaths that are not fw.
    for (set<int>::iterator j = this_walk.begin( ); j != this_walk.end( ); )
    {
        int u = *j;
        if ( !_local_unipaths_fw[u] )
            this_walk.erase(j++);
        else
            ++j;
    }

    if ( this_walk.empty() ) continue;
    n_walked++;
    walked[ ids[i] ] = True;
    
    // Mark these unipaths as walked.
    set<int>::iterator iter;
    //vec<int> walk_uIDs;
    //vecKmerPath walk_unipaths;
    for ( iter = this_walk.begin(); iter != this_walk.end(); iter++ ) {
      //walk_unipaths.push_back( _local_unipaths[*iter] );
      //walk_uIDs.push_back( *iter );
      unipaths_walked[*iter] = True;
    }
    
    // Mark these unipaths' adjacencies as traversed.
    set<int>::iterator iter2;
    for ( iter = this_walk.begin(); iter != this_walk.end(); iter++ )
      for ( iter2 = this_walk.begin(); iter2 != this_walk.end(); iter2++ ) {
	int u1 = *iter, u2 = *iter2;
	if ( _AG.HasEdge( u1, u2 ) ) {
	  // Don't mark an adjacency if it's already been marked,
	  // either forward or reverse (i.e., with RC'ed unipaths.)
	  if ( unipath_adjs.HasEdge( u1, u2 ) ) continue;
	  if ( unipath_adjs.HasEdge( to_rc[u2], to_rc[u1] ) ) continue;
	  
	  unipath_adjs.AddEdge( u1, u2 );
	}
      }
    
  }
  
  // Find all of the adjacencies implied by the unipath walks.
  //vecKmerPath walk_unipaths;
  //for ( int i = 0; i < _n_local_unipaths; i++ )
  //if ( unipaths_walked[i] )
  //  walk_unipaths.push_back( _local_unipaths[i] );
  
  // Using the adjacency digraph we've created, make the HKP of this walk.
  HyperKmerPath hkp;
  BuildUnipathAdjacencyHyperKmerPath( _K, unipath_adjs, _local_unipaths, hkp );
  
  
  // Clear the non-walked unipaths out of the HyperKmerPath.
  vec<int> edges_to_delete;
  for ( int i = 0; i < _n_local_unipaths; i++ )
    if ( !unipaths_walked[i] )
      edges_to_delete.push_back( i );
  hkp.DeleteEdges( edges_to_delete );
  hkp.RemoveDeadEdgeObjects( );
  hkp.RemoveEdgelessVertices( );

  // Find simple loops that contain copy number one stuff and flatten them.

  vec<simple_loop> loops;
  GetSimpleLoops( hkp, loops );
  for ( int j = 0; j < loops.isize( ); j++ )
  {    const simple_loop& S = loops[j];
       const KmerPath& L = hkp.EdgeObject( S.vv );
       int min_CN = 1000000000;
       for ( int r = 0; r < L.NSegments( ); r++ )
       {    const KmerPathInterval& I = L.Segment(r);
            vec<longlong> con;
            Contains( _pathsdb, I, con );
            for ( int s = 0; s < con.isize( ); s++ )
            {    const tagged_rpint& t = _pathsdb[ con[s] ];
                 if ( t.ReadId( ) < _n_reads ) continue;
                 int u = _cloud_unipaths[ t.ReadId( ) - _n_reads ].Uid( );
                 min_CN = Min( min_CN, (*_global_predicted_CNs)[u] );    }    }
       if ( min_CN > _ploidy ) continue;
       int z = hkp.N( );
       hkp.AddVertices(1);
       hkp.GiveEdgeNewToVx( S.vv, S.v, z );
       hkp.GiveEdgeNewFromVx( S.vw, S.v, z );    }
  
  // Dump selected pairs.

  if ( Member( _LOCAL_DUMP, String("SELECTED_PAIRS") ) )
  {    Ofstream( out, _outhead + ".selected_pairs.fasta" );
       for ( size_t i = 0; i < selected_pair_IDs.size( ); i++ )
       {    longlong pID = selected_pair_IDs[i];
            int sep = _global_pairs->sep(pID), sd = _global_pairs->sd(pID);
            int64_t global_id1 = _global_pairs->ID1(pID);
            int64_t global_id2 = _global_pairs->ID2(pID);
            int64_t local_id1 = BinPosition( _local_to_global, global_id1 );
            int64_t local_id2 = BinPosition( _local_to_global, global_id2 );
            const basevector &rd1 = _reads[local_id1], &rd2 = _reads[local_id2];
            rd1.Print( out, "select_" + ToString(i) + "_1 local_id=" 
                 + ToString(local_id1) + " global_id=" + ToString(global_id1)
                 + " sep=" + ToString(sep) + " dev=" + ToString(sd)
                 + " walked=" + ( walked[i] ? "yes" : "no" ) 
		 + ( redundant[i] ? " redundant" : "" ) );
            rd2.Print( out, "select_" + ToString(i) + "_2 local_id=" 
                 + ToString(local_id2) + " global_id=" + ToString(global_id2)
                 + " sep=" + ToString(sep) + " dev=" + ToString(sd)
                 + " walked=" + ( walked[i] ? "yes" : "no" ) 
		 + ( redundant[i] ? " redundant" : "" ) );    }    }
  
  return hkp;
}








/**
   Function: MergeInsertWalksLG

   Takes <HyperKmerPaths> obtained from walking long-insert pairs in
   one neighborhood, and merges them into a single HyperKmerPath for the
   neighborhood.  Adapted from MergeNeighborhoods for the LG pipeline.

   Inputs:

      hypers - HyperKmerPaths, one for each long-insert pair, representing
         the possible sequences in the middle of that pair

   Output: a HyperBasevector representing the possible
         sequences of the entire neighborhood.
*/
void
SeedNeighborhood::MergeInsertWalks( const vec<HyperKmerPath>& HKPs,
				    const vecKmerPath & new_unipaths,
				    const vecbasevector & new_unibases )
{
  // Create a database and KmerBaseBroker for the new unipaths.
  vecKmerPath new_unipaths_rc = new_unipaths;
  for ( size_t i = 0; i < new_unipaths_rc.size(); i++ )
    new_unipaths_rc[i].Reverse();
  vec<tagged_rpint> new_unipathsdb;
  CreateDatabase( new_unipaths, new_unipaths_rc, new_unipathsdb );
  KmerBaseBroker new_kbb;
  new_kbb.Initialize( _K, new_unibases, new_unipaths, new_unipaths_rc, new_unipathsdb );
  
  HyperKmerPath HKP_merged( _K, HKPs );
  if (_DUMP_UNMERGED)
  {    HyperBasevector hb( HKP_merged, new_kbb );
       Ofstream( out, _outhead + ".unmerged.fasta" );
       out << hb;
       Ofstream( outd, _outhead + ".unmerged.dot" );
       HKP_merged.PrintSummaryDOT0w( outd, True, False, True, 0, True );    }

  // Merge HyperKmerPaths using InternalMerge.
  // We attempt two merges.
  NegativeGapValidator ngv(&new_kbb);

  // First merge.
  InternalMerge( HKP_merged, &ngv, 5000, 4000 );
  InternalMerge( HKP_merged, &ngv, 5000, 3000 );
  InternalMerge( HKP_merged, &ngv, 5000, 2000 );
  InternalMerge( HKP_merged, &ngv, 5000, 1000 );
  InternalMerge( HKP_merged, &ngv, 5000, 500 );
  HKP_merged.Zipper( );

  // Second merge.
  InternalMerge( HKP_merged, &ngv, 3500, 700 );
  HKP_merged.Zipper( );

  // Attempt to glue replicated edge sequence together
  if (_NHOOD_GLUEPERFECT_SIZE > 0) {
    GluePerfects(HKP_merged, new_kbb, _NHOOD_GLUEPERFECT_SIZE);
  }

  // Housecleaning
  HKP_merged.ReduceLoops( );
  HKP_merged.CompressEdgeObjects( );
  HKP_merged.RemoveDeadEdgeObjects( ); 
  HKP_merged.RemoveEdgelessVertices( );    

  // Convert the HyperKmerPath to a HyperBasevector.
  _hbv = HyperBasevector( HKP_merged, new_kbb );
}





