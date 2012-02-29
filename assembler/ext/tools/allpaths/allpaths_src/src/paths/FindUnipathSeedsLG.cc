///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <pthread.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "Set.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/FindUnipathSeedsLG.h"
#include "paths/SeedNeighborhood.h"
#include "paths/simulation/Placement.h"

// MakeDepend: library PTHREAD

// A whole bunch of variables that are shared by parallel threads.

static const digraphE<fsepdev>* FG_ptr;
static const vec<int>* predicted_copyno_ptr;
static const vec<Bool>* branches_ptr;
static vec< vec<ustart> >* cloud_unipaths_ptr;
static vec<SeedNeighborhood *> * nhoods_ptr;

void *ThreadedFunction(void *threadid) {

  int *tid = (int *) threadid;
  int i = *tid;
  (*nhoods_ptr)[i]->FindUnipathCloud( (*FG_ptr), (*predicted_copyno_ptr), (*branches_ptr), 2000 );
  (*nhoods_ptr)[i]->ExpandUnipathCloud( (*FG_ptr), (*predicted_copyno_ptr), 2000 );
  (*cloud_unipaths_ptr)[i] = (*nhoods_ptr)[i]->GetCloudUnipaths();
  return 0;
}

template digraphE<fsepdev>::digraphE(digraphE<fsepdev> const&, vec<int> const&);






// These variables are only for logging purposes (via LogSeedStatus)
static vec<SeedStatus> * seed_status_ptr;
static ostream * seed_log_ptr;
static int max_uID_size;

// Helper function to maintain a log of all changes to seed statuses
// (i.e., unipaths being marked as non-seeds for any number of reasons.)
void LogSeedStatus( const int uID, const SeedStatus status, const String & reason ) {
  (*seed_status_ptr)[uID] = status;
  
  // The log contains 3 tab-delimited pieces of information:
  // 1. Unipath ID
  // 2. New seed status
  // 3. Reason for change
  String uID_str = ToString( uID );
  String status_str = status_names[status];
  // Pad the strings
  while ( uID_str.isize() < max_uID_size ) uID_str = "0" + uID_str;
  while ( status_str.size() < 23 ) status_str += " ";
  
  (*seed_log_ptr) << "UNIPATH#" << uID_str << "\t" << status_str << "\t" << reason << endl;
}





void FindUnipathSeeds(

     // inputs:

     const int K,
     const int ploidy,
     const int MIN_KMERS_IN_SEED,      // smallest seed to consider
     const int MIN_KMERS_IN_NEIGHBORS, // seed plus neighbors must be this big
     const vec<int>& ulen,             // unipath lengths
     const vec<unipath_id_t>& to_rc,   // unipath involutions
     const vec<int>& predicted_copyno, // predicted copy number for unipaths
     const digraphE<fsepdev>& FG,      // graph of all normal unipaths
     const Bool USE_TRUTH,             // if so, cheat by discarding RC unipaths
     VecPlacementVec locs,             // unipath locations on reference
     const int MAX_SEED_DIST,          // must use seed if farther than this
     const vec<Bool> & branches,
     const vecKmerPath* global_paths,
     const vecKmerPath* global_paths_rc,
     const vec<tagged_rpint>* global_pathsdb,
     const vec<int>* global_read_lengths,
     const vec<ReadLocationLG>* global_unilocs,
     const vec<longlong>* global_unilocs_index,
     const PairsManager* global_pairs,

     // output:

     vec<int>& seeds,
     vec<SeedStatus>& seed_status,
     ostream & seed_log,

     // extra params:

     const unsigned NUM_THREADS
          )
{
  int nuni = ulen.size( );
  seed_status.resize_and_set(nuni, SEED_GOOD);
  seed_status_ptr = &seed_status;
  seed_log_ptr = &seed_log;
  max_uID_size = ToString( nuni ).size();
  
  // Mark several types of unipaths as bad...
  for ( int v = 0; v < nuni; v++ ) {
    
    // Unipaths that are RC-on-ref (if a reference was given)
    if ( USE_TRUTH && !locs[v].empty( ) ) {
      Bool have_fw = False;
      for ( PlacementVec::size_type u = 0; u < locs[v].size( ); u++ ) {
	const placement& p = locs[v][u];
	if ( !p.Rc( ) ) have_fw = True;
      }
      if ( !have_fw )
	LogSeedStatus( v, SEED_RC_ON_REF, "[CHEATING] This unipath has no forward placements on the reference genome" );
    }
    
    // Unipaths with odd copy numbers
    if ( predicted_copyno[v] > ploidy )
      LogSeedStatus( v, SEED_HIGH_CN, "This unipath has a predicted CN of " + ToString( predicted_copyno[v] ) + " > ploidy" );
    else if ( predicted_copyno[v] == 0 )
      LogSeedStatus( v, SEED_ZERO_CN, "This unipath has a predicted CN of 0" );
    
    // Short unipaths
    if ( ulen[v] < MIN_KMERS_IN_SEED )
      LogSeedStatus( v, SEED_SHORT, "This unipath has a length of " + ToString( ulen[v] ) + " kmers (MIN_KMERS_IN_SEED = " + ToString( MIN_KMERS_IN_SEED ) + ")" );
    
    // Unipaths without connections in link graph
    if ( FG.From(v).empty( ) && FG.To(v).empty( ) )
      LogSeedStatus( v, SEED_ISOLATED, "This unipath has no connections in the link graph" );
  }

  // Eliminate seeds that aren't connected to much.

  for ( int v = 0; v < nuni; v++ ) {
    if ( seed_status[v] != SEED_GOOD ) continue;
    int total = ulen[v];
    for ( int j = 0; j < FG.From(v).isize( ); j++ )
      total += ulen[ FG.From(v)[j] ];
    for ( int j = 0; j < FG.To(v).isize( ); j++ )
      total += ulen[ FG.To(v)[j] ];
    if ( total < MIN_KMERS_IN_NEIGHBORS )
      LogSeedStatus( v, SEED_ISOLATED, "This unipath, plus its " + ToString(  FG.From(v).isize() + FG.To(v).isize() ) + " immediate neighbor(s) in the link graph, have a total kmer length of " + ToString( total ) + " (MIN_KMERS_IN_NEIGHBORS = " + ToString( MIN_KMERS_IN_NEIGHBORS ) + ")" );
  }
  
  // REDUNDANT SEED REMOVAL
  // Start by choosing all eligible unipaths as seeds.
  // Then go through them in some clever order,
  // and eliminate anything whose elimination would
  // not create an unaccceptably large gap.
  // METHOD CHANGED, SEE BELOW!!
  vec<int> seed_elim_order;
  
  // Elimination order possibilities:
  
  // 1. Try to drop small unipaths first
  WhatPermutation( ulen, seed_elim_order, less<int>(), false );
  
  // 2. Try to drop repetitive unipaths first, refined by length
  {

    // We would like to start with longest unipaths with (copy-number) cn = ploidy first.
    // However we would like to include unipaths with cn < ploidy before
    // those with cn > ploidy
    // vec< pair<float,int> > copyno_and_len( ulen.size() );
    // vec<float> fcn_tmp( predicted_copyno.size() );
    // int maxcn = Max( predicted_copyno );
  
    // for ( int k = 0; k < predicted_copyno.isize(); k++ ){
    //   if ( predicted_copyno[k] == 0 )
    //        fcn_tmp[k] = maxcn +1; /// will move cn = 0 (predicted) to the end
    //  else if ( predicted_copyno[k] < ploidy )
	// cn closer to ploidy first
    //    fcn_tmp[k] = ploidy + ( 1.0 - (float)predicted_copyno[k] / (float)ploidy );
    //   else fcn_tmp[k] = predicted_copyno[k];
    //  }
    // for(int i=0; i<nuni; i++)
    //   copyno_and_len[i] = make_pair( fcn_tmp[i] , -ulen[i] );
    // WhatPermutation( copyno_and_len, seed_elim_order, less< pair<float,int> >(), false );
  }

  cout << Date( ) << ": find redundant seeds by cloud building, <= 100 dots" << endl;
  vec< vec<int> > supportx(nuni), support(nuni);

  longlong count = 0, pass = 0, ngoods = 0;
  for( vec<int>::reverse_iterator uni_iter = seed_elim_order.rbegin();
       uni_iter != seed_elim_order.rend(); uni_iter++ ) {
    int v = *uni_iter;
    if ( seed_status[v] == SEED_GOOD ) ngoods += ulen[v];
  }
  
  vec<int> slim;
  for( vec<int>::reverse_iterator uni_iter = seed_elim_order.rbegin();
       uni_iter != seed_elim_order.rend(); uni_iter++ ) {
    slim.push_back( *uni_iter );
  }

  // The idea of what follows is to build a small neighborhood around each seed
  // and mark as redundant everything in the neighborhood except the seed.

  int mcn_other = ploidy;
  int NHOOD_RADIUS = 20000;

  int batch_size = NUM_THREADS;

  vec<int> batch;
  vec< vec<ustart> > cloud_unipaths;

  FG_ptr = &FG;
  predicted_copyno_ptr = &predicted_copyno;
  branches_ptr = &branches;
  cloud_unipaths_ptr = &cloud_unipaths;

  for ( int si = 0; si < slim.isize( ); si++ ) {
    batch.clear( );
    for ( ; si < slim.isize( ); si++ ) {
      int v = slim[si];
      if ( seed_status[v] != SEED_GOOD ) continue;
      for ( int j = 0; j < ulen[v]; j++ )
	if ( count++ % (ngoods/100) == 0 ) Dot( cout, pass++ );
      batch.push_back(v);
      if ( batch.isize( ) == batch_size ) break;
    }
    si--;

    cloud_unipaths.clear_and_resize( batch.size( ) );

    int JOBS = batch.size( );
    // Create a set of SeedNeighborhood objects, one for each seed.
    vec<SeedNeighborhood *> nhoods;
    for ( int i = 0; i < JOBS; i++ ) {
      nhoods.push_back( new SeedNeighborhood(
					     K, ploidy, mcn_other, global_paths,
					     global_paths_rc, global_pathsdb, &ulen, global_read_lengths, &predicted_copyno,
					     global_pairs, global_unilocs, global_unilocs_index ) );
      nhoods.back()->SetSeedID( batch[i] );
    }
    SeedNeighborhood::_NHOOD_RADIUS          = NHOOD_RADIUS;
    SeedNeighborhood::_NHOOD_RADIUS_INTERNAL = NHOOD_RADIUS;
    nhoods_ptr = &nhoods;
    
    // Process the SeedNeighborhoods in parallel, creating each seed's
    // unipath cloud via the ThreadedFunction call.
    vec<pthread_t> threads(JOBS);
    vec<unsigned int> ids(JOBS);
    for ( int i = 0; i < JOBS; i++ ) {
      ids[i] = i;
      pthread_create( &threads[i], NULL,
		      ThreadedFunction, (void*) &ids[i] );
    }
    for ( int i = 0; i < JOBS; i++ )
      pthread_join( threads[i], NULL );
	  
    for ( int i = 0; i < JOBS; i++ )
      delete nhoods[i];

    for ( int i = 0; i < batch.isize( ); i++ ) {
      int v = batch[i];
      if ( seed_status[v] != SEED_GOOD ) continue;
      for ( int j = 0; j < cloud_unipaths[i].isize( ); j++ ) {
	int u = cloud_unipaths[i][j].Uid( );
	if ( Abs( cloud_unipaths[i][j].Start( ) ) > MAX_SEED_DIST )
	  continue;
	supportx[u].push_back(v), support[v].push_back(u);
	if ( u != v ) {
	  if ( seed_status[u] == SEED_GOOD ) {
	    ngoods -= ulen[u];
	    count -= ulen[u];
	    int actual_pass = (count-1) / (ngoods/100) + 1;
	    for ( ; pass < actual_pass; pass++ )
	      Dot( cout, pass++ );
	  }
	  ostringstream reason;
	  reason << "This unipath appears in the neighborhood of good unipath " << v << " at a distance of " << Abs( cloud_unipaths[i][j].Start() ) << " (MAX_SEED_DIST = " << MAX_SEED_DIST << ") and is thus redundant";
	  LogSeedStatus( u, SEED_REDUNDANT, reason.str() );
	}
      }
    }
  }
 
  // Now remove seeds if they do not add to support.

  cout << "\n" << Date( ) << ": removing seeds that do not add to support" << endl;
  for ( int v = 0; v < nuni; v++ )
    UniqueSort( supportx[v] );
  for ( int v = 0; v < nuni; v++ ) {
    if ( seed_status[v] != SEED_GOOD ) continue;
    
    // Check whether there are any other seeds whose unipath neighborhood is
    // a superset of this seed's unipath neighborhood.
    vec<int> common = supportx[ support[v][0] ];
    for ( int l = 1; l < support[v].isize( ); l++ )
      common = Intersection( common, supportx[ support[v][l] ] );
    if ( common.size( ) > 1 ) {
      
      String reason = "This unipath's neighborhood is a subset of the following other unipaths' neighborhoods:";
      for ( int i = 0; i < common.isize(); i++ )
	reason += " " + ToString( common[i] );
      LogSeedStatus( v, SEED_REDUNDANT, reason );
      
      for ( int l = 0; l < support[v].isize( ); l++ )
	supportx[ support[v][l] ].EraseValue(v);
    }
  }
  
  // If a unipath and its RC are currently both chosen as seeds,
  // remove the higher-numbered unipath from the list.

  for ( int v = 0; v < nuni; v++ )
    if ( seed_status[v]          == SEED_GOOD &&
	 seed_status[ to_rc[v] ] == SEED_GOOD &&
	 v > to_rc[v] )
      LogSeedStatus( v, SEED_RC_OF_SEED, "This unipath is the RC of unipath #" + ToString( to_rc[v] ) + ", which is being chosen as a seed" );
  
  
  // Now record which seeds were chosen

  seeds.clear( );
  for(int v=0; v < nuni; v++)
    if( seed_status[v] == SEED_GOOD ) {
      LogSeedStatus( v, SEED_GOOD, "Seed #" + ToString( seeds.isize() ) );
      seeds.push_back(v);
    }
  
  // Report on the reasons for rejecting unipaths from consideration as seeds
  cout << "Reasons for rejection as potential unipaths..." << endl;
  cout << "There are a total of " << nuni << " unipaths." << endl;
  cout << "Rejections are applied in the following order:" << endl;
  cout << "Isolated in link graph:\t"
       << seed_status.CountValue( SEED_ISOLATED ) << endl;
  cout << "Shorter than MIN_KMERS:\t"
       << seed_status.CountValue( SEED_SHORT ) << endl;
  cout << "Copy number > ploidy:\t"
       << seed_status.CountValue( SEED_HIGH_CN ) << endl;
  cout << "Copy number is zero:\t"
       << seed_status.CountValue( SEED_ZERO_CN ) << endl;
  cout << "RC on ref (cheat!):\t"
       << seed_status.CountValue( SEED_RC_ON_REF ) << endl;
  cout << "Deemed redundant:\t"
       << seed_status.CountValue( SEED_REDUNDANT ) << endl;
  cout << "RC of another seed:\t"
       << seed_status.CountValue( SEED_RC_OF_SEED ) << endl;
  cout << "Good!\t\t\t" << seed_status.CountValue( SEED_GOOD ) << endl;
  
  cout << endl << "Ultimately selected " << seeds.size( ) << " seeds." << endl;
}








// Evaluate the seed selection process by aligning the seed unipaths to
// reference and measuring the gaps between them
void
EvalUnipathSeeds( const vec<int> & seeds, const vec<int>& ulen,
		  const vecbasevector & genome, const VecPlacementVec& locs,
		  const int MAX_SEED_DIST, const vec<int> & predicted_copyno,
		  const digraphE<fsepdev> & FG, const vec<SeedStatus>& seed_status ) {

  int nuni = ulen.size( );
  
  vec<int> seed_lengths( seeds.isize( ) );
  for ( int i = 0; i < seeds.isize( ); i++ )
    seed_lengths[i] = ulen[ seeds[i] ];
  
  // Cancel out if no truth is supplied
  if ( !genome.size( ) || !locs.size( ) ) return;
  
  
  
  // Tabulate the true locations of seeds.
  vec< pair<placement,int> > P;
  for ( int i = 0; i < seeds.isize( ); i++ ) {
    int v = seeds[i];
    for ( PlacementVec::size_type j = 0; j < locs[v].size( ); j++ )
      P.push_back( make_pair( locs[v][j], v ) );
  }
  Sort(P);
  
  vec< vec<ho_interval> > cov( genome.size( ) );
  for ( int i = 0; i < P.isize( ); i++ )
    cov[ P[i].first.GenomeId( ) ].push_back(ho_interval( P[i].first.pos( ), P[i].first.Pos( ) ) );
  
  
  vec<ho_interval> un;
  vec<int> un_g;
  vec<int> un_length, un_gap_length, un_end_length, un_contig_length;
  // Find covered and uncovered stretches of the reference contigs
  for ( size_t g = 0; g < genome.size( ); g++ ) {
    vec<ho_interval> ung;
    Uncovered( genome[g].size( ), cov[g], ung );
    un.append(ung);
    
    for ( size_t j = 0; j < ung.size( ); j++ ) {
      un_g.push_back(g);
      int len = ung[j].Length( );
      un_length.push_back( len );
      if ( j == 0 && ung.size( ) == 1 )
	un_contig_length.push_back( len );
      // If this stretch is at either end of the reference contig,
      // it's a 'missing end'; otherwise it's a 'gap'
      else if ( j == 0 || j + 1 == ung.size( ) )
	un_end_length.push_back( len );
      else
	un_gap_length.push_back( len );
    }
  }
  ReverseSortSync( un_length, un_g, un );
  vec< vec<String> > rows;
  for ( size_t j = 0; j < Min( un.size( ), 20UL ); j++ ) {
    vec<String> row;
    row.push_back( ToString( un_length[j] ) );
    row.push_back( ToString( un_g[j] ) + "." + ToString( un[j].Start( ) )
		   + "-" + ToString( un[j].Stop( ) ) );
    rows.push_back(row);
  }
  
  cout << "\nLargest uncovered stretches:\n";
  PrintTabular( cout, rows, 2, "rl" );
  flush(cout);
  
  
  // Print out N50 sizes of seeds, gaps, and missing ends
  Sort( seed_lengths );
  Sort( un_gap_length );
  Sort( un_end_length );
  Sort( un_contig_length );
  cout << endl;
  cout << "N50 sizes (in base space):" << endl;
  cout << "Seed unipaths:                       ";
  if ( seed_lengths.nonempty( ) ) cout << N50( seed_lengths );
  else cout << "N/A";
  cout << endl;
  if ( un_gap_length.isize( ) )
    cout << "Gaps between seed unipaths:          "
	 << N50( un_gap_length ) << endl;
  else
    cout << "There are no gaps between seed unipaths." << endl;
  if ( un_end_length.isize( ) )
    cout << "Uncovered ends on reference contigs: "
	 << N50( un_end_length ) << endl;
  else
    cout << "There are no uncovered ends on reference contigs." << endl;
  if ( un_contig_length.isize( ) )
    cout << "Totally uncovered reference contigs: "
	 << N50( un_contig_length ) << endl;
  else
    cout << "There are no totally uncovered reference contigs." << endl;
  cout << endl;
  
  // Print seeds and statuses in order along the genome, if there aren't
  // too many of them.
  if ( seeds.size( ) > 100 ) return;
  
  cout << "\nWalking along the reference:" << endl;
  for ( int p=0; p<P.isize(); p++ )
    if ( seed_status[P[p].second] == SEED_GOOD )
      cout << P[p].first << '\t'
	   << P[p].second << '\t'
	   << locs[ P[p].second ].size() << " locs\t"
	   << status_names[seed_status[P[p].second]] << '\n';
  
  cout << endl;
}
