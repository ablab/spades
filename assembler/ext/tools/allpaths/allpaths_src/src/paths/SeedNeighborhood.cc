///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Set.h"
#include "VecUtilities.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "paths/ReadsToPathsCoreX.h" // ReadsToPathsCoreY
#include "paths/SeedNeighborhood.h"
#include "paths/Unipath.h"
#include "feudal/Algorithms.h"
#include "TaskTimer.h"



int SeedNeighborhood::_NHOOD_RADIUS;
int SeedNeighborhood::_NHOOD_RADIUS_INTERNAL;
int SeedNeighborhood::_NHOOD_GLUEPERFECT_SIZE;
int SeedNeighborhood::_NHOOD_VERBOSITY;
Bool SeedNeighborhood::_USE_ACYCLIC;
Bool SeedNeighborhood::_DUMP_UNMERGED;
Bool SeedNeighborhood::_DUMP_LOCAL_UNIBASES;

/* Constructor: immediately loads all global data structures
 *
 * REQUIRED_INPUT:
 * seed_unipath_ID: The ID of the seed unipath (in global unipath numbering)
 * K: Kmer size that has been used for global pathing
 * paths: All reads, ordered by global read ID, pathed together
 * unipath_lengths: Lengths of global unipaths
 * pairs: All read pairs; reads are given global read IDs
 * pairs_index: Index of (global) read ID to (global) pair ID
 * unilocs: Locations of (global) reads on (global) unipaths
 * unilocs_index: Index of (global) unipaths to (global) read IDs
 *
 ******************************************************************************/
SeedNeighborhood::SeedNeighborhood
( const int K,
  const int ploidy,
  const int mcn_other,
  const vecKmerPath * paths,
  const vecKmerPath * paths_rc,
  const vec<tagged_rpint> * pathsdb,
  const vec<int> * unipath_lengths,
  const vec<int> * read_lengths,
  const vec<int> * predicted_CNs,
  const PairsManager * pairs,
  const vec<ReadLocationLG> * unilocs,
  const vec<longlong> * unilocs_index,
  const vec<String>& LOCAL_DUMP )
  : _K( K ),
    _ploidy( ploidy ),
    _mcn_other( mcn_other ),
    _global_paths( paths ),
    _global_paths_rc( paths_rc ),
    _global_pathsdb( pathsdb),
    _global_unipath_lengths( unipath_lengths ),
    _global_read_lengths ( read_lengths ),
    _global_predicted_CNs( predicted_CNs),
    _global_pairs( pairs ),
    _global_unilocs( unilocs ),
    _global_unilocs_index( unilocs_index ),
    _LOCAL_DUMP( LOCAL_DUMP )
{
  _seed_ID = -1;
  
  _n_cloud_unipaths = _n_reads = _n_pairs = _n_local_unipaths = 0;
  
  // All logging output is ignored until a call is made to
  // SeedNeighborhood::SetLog( ) later
  _log = new ofstream( "/dev/null" );

  // Similarly, eval_subdir is only created upon request (by
  // this->SetEvalSubdir).
  _eval_subdir = NULL;
  
  // Set counters to 0.
  _n_inserts_walked = _n_inserts_total = 0;
  _T_pathing  = _T_insert_walking  = _T_insert_merging = 0;
}



SeedNeighborhood::~SeedNeighborhood( )
{
  if ( _log != &cout ) { delete _log; _log = NULL; }
  if ( _eval_subdir ) { delete _eval_subdir; _eval_subdir = NULL; }
}



// Optional: Set up output file (or stdout) for logging
void
SeedNeighborhood::SetLog( const String & logfile )
{
  if ( _log != &cout ) { delete _log; _log = NULL; }
  
  if ( logfile == "stdout" )
    _log = &cout;
  else if ( logfile == "" )
    _log = new ofstream( "/dev/null" );
  else
    _log = new ofstream( logfile.c_str( ) );
  
  LogDate( "Begin logging output for seed ID " + ToString( _seed_ID ), 2 );
}

// Optional: Create an eval_subdir subdirectory, that can be used to
// dump possibly large intermediary files.
void SeedNeighborhood::SetEvalSubdir( const String &evaldir )
{
  _eval_subdir = new String( evaldir );
  Mkpath( *_eval_subdir );
}


/* Make the cloud of unipaths in the vicinity of this seed
 * Note that these "cloud unipaths" are entirely different from the "local
 * unipaths" which we will create later
 *
 * INPUT: Unipath link graphs (made by BuildUnipathLinkGraphs) and copy numbers
 * MUST ALREADY BE RUN: Only the constructor
 *
 * METHOD:
 * 1. Run FindUnipathCloud to find unipaths in the neighborhood of this seed.
 * 2. Run ExpandUnipathCloud to add some unipaths with predicted copy number
 *    <= ploidy.
 * 
 ******************************************************************************/
void
SeedNeighborhood::MakeUnipathCloud
( const digraphE<fsepdev> & LG, const vec<int> & predicted_CNs,
  const vec<Bool> & branches )
{
  LogDate( "Extending seed's neighborhood to nearby unipaths..." );
  
  // Find the unipaths in the neighborhood of this seed
  const int MAX_DEV = 2000;
  FindUnipathCloud( LG, predicted_CNs, branches, MAX_DEV );
  // Add in unipaths of predicted copy number <= ploidy
  ExpandUnipathCloud( LG, predicted_CNs, MAX_DEV );

  _n_cloud_unipaths = _cloud_unipaths.isize( );
  LogDate( "Neighborhood now contains " + ToString( _n_cloud_unipaths )
	   + " unipaths." );
  
  // Save ids of cloud unipaths (sorted).
  if ( _eval_subdir ) {
    vec<int> ids;
    ids.reserve( _n_cloud_unipaths );
    for (int ii=0; ii<_n_cloud_unipaths; ii++)
      ids.push_back( _cloud_unipaths[ii].Uid( ) );
    sort( ids.begin( ), ids.end( ) );
    
    String ids_file = *_eval_subdir + "/cloud_unipaths.ids";
    ofstream out( ids_file.c_str( ) );
    for (int ii=0; ii<ids.isize( ); ii++)
      out << ids[ii] << "\n";
    out.close( );
  }
}


/* Find the reads and pairs in this neighborhood
 *
 * INPUT: Bases in reads
 * MUST ALREADY BE RUN: MakeUnipathCloud
 *
 * METHOD:
 * 1. Run FindPrimaryReadCloud to find all reads (and unibases) that
 *    align to these unipaths
 * 2. Create a vecbasevector of local reads
 * 
 ******************************************************************************/
void
SeedNeighborhood::MakeReadCloud( const vecbvec &global_reads_bases,
				 const vecbvec &unibases, const Bool& LOCAL_PRIMARY )
{
  LogDate( "Populating neighborhood with reads..." );
  ForceAssert( _n_cloud_unipaths != 0 );
  
  FindPrimaryReadCloud( global_reads_bases, LOCAL_PRIMARY );

  size_t primary_count = _read_IDs.isize();
  LogDate( "Primary read cloud contains " + ToString( primary_count )
	   + " reads." );
  
  FindSecondaryReadCloud( );

  size_t secondary_count = _read_IDs.isize() - primary_count;
  LogDate( "Secondary read cloud contains " + ToString( secondary_count )
	   + " reads." );

  
  // Save ids of cloud reads (both primary and secondary, sorted).

  if ( _eval_subdir ) {
    String ids_file = *_eval_subdir + "/cloud_reads.ids";
    ofstream out( ids_file.c_str( ) );
    for (int ii=0; ii<_read_IDs.isize( ); ii++)
      out << _read_IDs[ii].first << "\n";
    out.close( );
  }

  // Basic report on the reads in this cloud

  _n_reads = _read_IDs.isize( );
  _n_pairs = _pair_IDs.isize( );
  LogDate( "Neighborhood contains " + ToString( _n_reads )
	   + " reads and " + ToString( _n_pairs ) + " pairs." );

  // Array to map read IDs in this neighborhood to global read IDs
  _local_to_global.resize( _n_reads, -1 );
  for ( int i = 0; i < _n_reads; i++ )
    _local_to_global[ i ] = _read_IDs[i].first;

  LogDate( "Global/Local read ID lookup tables computed." );

  // Reserve memory for local reads
  int n_objects = _n_reads + _n_cloud_unipaths;
  longlong n_bases = 0;
  for ( int i = 0; i < _n_reads; i++ )
    n_bases += (*_global_read_lengths)[ _local_to_global[i] ];
  for (int ii=0; ii<_n_cloud_unipaths; ii++)
    n_bases += unibases[ _cloud_unipaths[ii].Uid( ) ].size( );

  LogDate( "Reserving space for local reads." );
  
  _reads.clear( );
  _reads.Reserve( n_bases / 16 + n_objects , n_objects );

  LogDate( "Extracting local reads." );
  
  // Create the vecbasevector of local reads.
  int n_frag_reads = 0;
  _reads_fw.clear( );
  for ( int i = 0; i < _n_reads; i++ ) {
    basevector read = global_reads_bases[ _local_to_global[i] ];

    _reads_fw.push_back(read);
    if ( _read_IDs[i].second == ORIENT_RC ) _reads_fw.back( ).ReverseComplement( );
    
    if ( _global_pairs->isUnpaired( _local_to_global[i] ) ) {
      n_frag_reads++;
      // Flip unpaired reads so they face forward in the local neighborhood.
      if ( _read_IDs[i].second == ORIENT_RC ) {
	read = read.ReverseComplement();
	_read_IDs[i].second = ORIENT_FW;
      }
    }
      
    _reads.push_back( read );
  }

  LogDate( "Adding local unibases to reads." );
  
  // Add unibases. WARNING! This is done as a last step, so they appear
  // as the last objects in <reads>.  Notice also that the local_to_global map
  // is not defined here.
  for (int ii=0; ii<_n_cloud_unipaths; ii++)
    _reads.push_back( unibases[ _cloud_unipaths[ii].Uid( ) ] );
  for (int ii=0; ii<_n_cloud_unipaths; ii++)
    _reads_fw.push_back( unibases[ _cloud_unipaths[ii].Uid( ) ] );

  LogDate( ToString( _n_reads ) + " reads, "
	   + ToString( _n_cloud_unipaths ) + " unipaths, and "
	   + ToString( _n_pairs ) + " pairs into memory." );
	   
  // Report on the number of [fragment/jumping] reads [with/without] partner
  // in the local assembly.
  int n_reads_no_partner = _n_reads - 2*_n_pairs;
  
  LogDate( "[Paired reads: " + ToString( _n_pairs * 2 ) + " w/partner, " +
	   ToString( n_reads_no_partner - n_frag_reads) + " w/o]" );
  LogDate( "[Fragment reads: " + ToString( n_frag_reads ) + "]" );
  LogDate( "[Unibases: " + ToString( _n_cloud_unipaths ) + "]" );
}


/* Assemble this seed's local reads into a unipath graph and HyperKmerPath
 * Note that the "local unipaths" created here are entirely different from the
 * "cloud unipaths" found in the vicinity of the seed unipath
 *
 * INPUT: None
 * MUST ALREADY BE RUN: MakeReadCloud
 *
 * METHOD:
 * 1. Path the reads - this creates a new local kmer numbering system
 * 2. Create a KmerBaseBroker for the new kmers
 * 3. Create unipaths out of the paths
 * 4. Create a HyperKmerPath out of the unipaths
 *
 ******************************************************************************/
void
SeedNeighborhood::MakeLocalAssembly( const vecbasevector& unibases,
     const vec< vec<int> >& unibases_next, const vec<int>& to_rc, const int pass )
{
  LogDate( "Pathing " + ToString( _n_reads ) + " reads..." );
  if ( _n_reads == 0 ) {
    *_log << "\tThe local read cloud is empty.  No pathing to be done." << endl;
    return;
  }
  
  // Create read paths (this is a timed operation.)
  double t_start = WallClockTime();
  
  _paths.clear().shrink_to_fit();
  _paths_rc.clear().shrink_to_fit();
  _pathsdb.clear_and_resize(0);
  ReadsToPathsCoreY( _reads, _K, _paths, _paths_rc, _pathsdb );
  _T_pathing = WallClockTime() - t_start;
  

  LogDate( "Assembling neighborhood locally..." );
  
  // Create KmerBaseBroker for this neighborhood

  _kbb.Initialize( _K, _reads, _paths, _paths_rc, _pathsdb );
  
  // Create reconstructed unipaths
  // NOTE: These are *not* the same as the cloud unipaths that were originally
  // found in the vicinity of the seed!
  _local_unipaths.clear().shrink_to_fit();
  _local_unipathsdb.clear_and_resize(0);
  Unipath( _paths, _paths_rc, _pathsdb, _local_unipaths, _local_unipathsdb );
  _n_local_unipaths = _local_unipaths.size( );
  /*
  PRINT(pass); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  PRINT(_n_local_unipaths); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  PRINT( _local_unipathsdb.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  */
  
  LogDate( "Reads reconstructed into " + ToString( _n_local_unipaths ) + " unipaths." );

  // Determine which local unipaths are in _reads_fw.  Very bad.

  vecbvec kmers;
  kmers.reserve(kmerCount(_reads_fw.begin(),_reads_fw.end(),_K));

  basevector b;
  typedef vecbvec::iterator Itr;
  for ( Itr itr(_reads_fw.begin()), end(_reads_fw.end()); itr != end; ++itr )
  {
      bvec const& bv = *itr;
      bvec::size_type nnn = bv.size();
      if ( nnn >= static_cast<unsigned>(_K) )
      {
          nnn -= _K;
          for ( bvec::size_type j = 0; j <= nnn; ++j )
              kmers.push_back(b.SetToSubOf(bv,j,_K));
      }
  }

  kmers.UniqueSort();

  // Figure out which unipaths are in the neighborhood, as opposed to being
  // in either the neighborhood or its reverse complement.  Somehow, in an
  // incarnation of the algorithm from years ago in LocalizeReads.cc, only the 
  // 'forward' stuff was present.
  _local_unipaths_fw.resize_and_set( _local_unipaths.size( ), False );
  for ( uint i = 0; i < _local_unipaths.size( ); i++ )
  {    basevector u = _kbb.Seq( _local_unipaths[i] );
       int present = 0, total = 0;
       for ( int j = 0; j <= u.isize( ) - _K; j++ )
       {    total++;
            b.SetToSubOf( u, j, _K );
            if ( BinPosition(kmers,b) != -1 ) present++;    }
       if ( present == total ) _local_unipaths_fw[i] = True;    }
  
  // Build local HyperKmerPath

  BuildUnipathAdjacencyGraph( _paths, _paths_rc, _pathsdb, _local_unipaths, _local_unipathsdb, _AG );
  BuildUnipathAdjacencyHyperKmerPath( _K, _AG, _local_unipaths, _hkp );

  ForceAssert( _hkp.EdgeObjectCount( ) == _n_local_unipaths ); // sanity check
  ForceAssert( _hkp.Edges( ) == VecOfKmerPath( _local_unipaths ) );

  // NEW CODE.
  
  if ( pass == 2 ) return;

  this ->AddGlobalConnections( unibases, unibases_next, to_rc );
  
  LogDate( "Unipaths built into a local HyperKmerPath." );
}

/**
 * Experimental code to add global connections to _reads. This is
 * called by MakeLocalAssemly( ) only if pass == 1.
 */
void
SeedNeighborhood::AddGlobalConnections( const vecbasevector& unibases,
					const vec< vec<int> >& unibases_next,
					const vec<int>& to_rc )
{
     // First convert each cloud unipath into a sequence of one or more local
     // unipaths.

     // cout << "\nCLOUD TO LOCAL MAP\n";
     vec< vec<int> > cloud_to_local(_n_cloud_unipaths);
     vec<longlong> con;
     for ( int i = 0; i < _n_cloud_unipaths; i++ )
     {    size_t r = _n_reads + i;
          const KmerPath& p = _paths[r];
          for ( int j = 0; j < p.NSegments( ); j++ )
          {    const KmerPathInterval& I = p.Segment(j);
               Contains( _local_unipathsdb, I, con );
               for ( int l = 0; l < con.isize( ); l++ )
               {    const tagged_rpint& t = _local_unipathsdb[ con[l] ];
                    if ( t.Rc( ) ) continue;
                    if ( cloud_to_local[i].empty( )
                         || cloud_to_local[i].back( ) != t.PathId( ) )
                    {    cloud_to_local[i].push_back( t.PathId( ) );    }    }    }
          // cout << "cloud unipath " << i << "["  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          //      << _cloud_unipaths[i].Uid( ) << "]:"; // XXXXXXXXXXXXXXXXXXXXXXXXX
          // for ( int j = 0; j < cloud_to_local[i].isize( ); j++ ) // XXXXXXXXXXXXX
          //      cout << " " << cloud_to_local[i][j]; // XXXXXXXXXXXXXXXXXXXXXXXXXX
          // cout << "\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               }

     // Now figure out how they are connected.

     // cout << "\nMAPS TO:\n";
     vec<int> to_left, to_right;
     _hkp.ToLeft(to_left), _hkp.ToRight(to_right);
     const int max_pile = 50;
     const int max_dist = 2000;
     vec<Bool> reaches_right( _n_cloud_unipaths, False );
     set<int> used;
     vec< pair<int,int> > X, Y;
     for ( int i = 0; i < _n_cloud_unipaths; i++ )
     {    if ( cloud_to_local[i].empty( ) ) continue;
          X.clear( ), Y.clear( );
          X.push( cloud_to_local[i].back( ), 0 );
          while( X.nonempty( ) && X.isize( ) + Y.isize( ) <= max_pile )
          {    pair<int,int> p = X.back( );
               X.pop_back( );
               Y.push_back(p);
               int v = to_right[ p.first ];
               for ( int j = 0; j < _hkp.From(v).isize( ); j++ )
               {    int e = _hkp.EdgeObjectIndexByIndexFrom( v, j );
                    int d = p.second + _hkp.EdgeObject(e).KmerCount( );
                    if ( d <= max_dist && !Member( used, e ) ) 
                    {    X.push( e, d );
                         used.insert(e);
                         for ( int j = 0; j < _n_cloud_unipaths; j++ )
                         {    if ( Member( cloud_to_local[j], e ) )
                              {    reaches_right[i] = True;
                                   // cout << i << " --> " << j << endl;
                                   break;    }    }    }
                    if ( reaches_right[i] ) break;    }
               if ( reaches_right[i] ) break;    }    }

     // Check the global unipaths.

     // cout << "\nGLOBAL CONNECTIONS\n";
     vec<int> extra_globals;
     used.clear( );
     vec<int> extra;
     vec<int> Z, W;
     for ( int i = 0; i < _n_cloud_unipaths; i++ )
     {    if ( reaches_right[i] ) continue;
          X.clear( ), Y.clear( ), Z.clear( ), W.clear( );
          X.push( _cloud_unipaths[i].Uid( ), 0 );
          Bool found = False;
          while( X.nonempty( ) && X.isize( ) + Y.isize( ) <= max_pile )
          {    pair<int,int> p = X.back( );
               X.pop_back( );
               Y.push_back(p);
               int u = p.first;
               for ( int j = 0; j < unibases_next[u].isize( ); j++ )
               {    int v = unibases_next[u][j];
                    int d = p.second + unibases[v].isize( ) - _K + 1;
                    if ( d <= max_dist && !Member( used,v ) )
                    {    X.push( v, d );
                         used.insert(v);
                         for ( int l = 0; l < _n_cloud_unipaths; l++ )
                         {    if ( v == _cloud_unipaths[l].Uid( ) )
                              {    found = True;
                                   // cout << "\n" << i << " --> " << l // XXXXXXXXX
                                   //      << endl; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                                   // PRINT2( X.size( ), Y.size( ) ); // XXXXXXXXXXX
                                   Z.push_back(v);
                                   vec<int> XY;
                                   for ( int m = 0; m < X.isize( ); m++ )
                                        XY.push_back( X[m].first );
                                   for ( int m = 0; m < Y.isize( ); m++ )
                                        XY.push_back( Y[m].first );
                                   UniqueSort(XY);
                                   while( Z.nonempty( ) )
                                   {    int v = Z.back( );
                                        Z.pop_back( );
                                        W.push_back(v);
                                        int rv = to_rc[v];
                                        for ( int m = 0; 
                                             m < unibases_next[rv].isize( ); m++ )
                                        {    int u = to_rc[ unibases_next[rv][m] ];
                                             if ( BinMember( XY, u )
                                                  && !Member( Z, u )
                                                  && !Member( W, u )
                                                  && u != _cloud_unipaths[i].Uid( ) )
                                             {    Z.push_back(u);    }    }    }
                                   W.push_back( _cloud_unipaths[i].Uid( ) );
                                   UniqueSort(W);
                                   break;    }    }
                         if (found) break;    }
                    if (found) break;    }
               if (found) break;    }
          if (found)
          {    
               /*
               extra.clear( );
               Y.append(X);
               for ( int j = 0; j < Y.isize( ); j++ )
                    extra.push_back( Y[j].first );
               UniqueSort(extra);
               extra_globals.append(extra);    
               */
               extra_globals.append(W);
                    }    }
     UniqueSort(extra_globals);
     // PRINT( extra_globals.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     // Push back the extra globals.

     int n_objects = _n_reads + _n_cloud_unipaths;
     longlong n_bases = 0;
     for ( int i = 0; i < _n_reads; i++ )
          n_bases += (*_global_read_lengths)[ _local_to_global[i] ];
     for (int ii=0; ii<_n_cloud_unipaths; ii++)
          n_bases += unibases[ _cloud_unipaths[ii].Uid( ) ].size( );
     for ( int i = 0; i < extra_globals.isize( ); i++ )
     {    int u = extra_globals[i];
          for ( int j = 0; j < unibases_next[u].isize( ); j++ )
          {    int v = unibases_next[u][j];
               n_objects++;
               n_bases += unibases[u].size( ) + 1;    }    }
     _reads.Reserve( n_bases / 16 + n_objects , n_objects );

     basevector bx;
     for ( int i = 0; i < extra_globals.isize( ); i++ )
     {    int u = extra_globals[i];
          for ( int j = 0; j < unibases_next[u].isize( ); j++ )
          {    int v = unibases_next[u][j];
               bx.resize( unibases[u].size( ) + 1 );
               for ( int l = 0; l < unibases[u].isize( ); l++ )
                    bx.Set( l, unibases[u][l] );
               bx.Set( unibases[u].size( ), unibases[v][_K-1] );
               _n_reads++;
               _reads.push_back(bx);    
               _reads_fw.push_back_reserve(bx);    }    } 

     // _n_reads = _reads.size( );
}

/* Walk inserts in this neighborhood with the help of fragments.  This converts
 * the local-assembly HyperKmerPath into a HyperBasevector, which represents
 * ONLY the parts of the HyperKmerPath containing successfully walked inserts.
 *
 * INPUT: Unibases corresponding to global unipaths
 * MUST ALREADY BE RUN: MakeLocalAssembly
 *
 * METHOD:
 * 1. Place all local reads and pairs onto the HyperKmerPath.
 * 2. Call WalkLongInserts to find pair closures.
 * 3. Call MergeInsertWalksLG to convert the set of pair closures to a single
 *    HyperBasevector, representing this neighborhood.
 * 4. Clean up the HyperBasevector with a call to ClipHKP.
 *
 * OUTPUT: returns 'false' on irregular exit (e.g. timeout)
 ******************************************************************************/
bool
SeedNeighborhood::MakeInsertWalks( )
{
  // Time out is commented out.
  TaskTimer timer;
  timer.StartWithTimeOut( 60 * 60 );  // time in seconds

  LogDate( "Aligning reads and pairs onto local HyperKmerPath..." );
  if ( _n_pairs == 0 ) {
    *_log << "\tThere are no paired reads in this neighborhood.  No insert walking to be done." << endl;
    return true;
  }
  ForceAssert( _n_reads != 0 );
  ForceAssert( _n_local_unipaths != 0 );
  
  const int min_size = 500; // minimum component size
 
  // Select a small subset of the pairs to walk using WalkLongInserts.
  // This is a heuristic to improve performance - because insert coverage is
  // higher than needed, and insert walking and merging are time-consuming and
  // can result in duplication.

  const int n_pairs_to_select = 100;
  int n_logical_pairs;
  vec<longlong> selected_pair_IDs = SelectPairsToWalk( n_pairs_to_select, n_logical_pairs );
  
  // The number of selected pairs may be less than n_pairs_to_select, if not
  // enough pairs in this neighborhood are logical.
  _n_inserts_total = selected_pair_IDs.isize( );
  
  LogDate( "SelectPairsToWalk has selected " + ToString( _n_inserts_total )
	   + " pairs out of " + ToString( _n_pairs )
	   + " total (" + ToString( n_logical_pairs ) + " logical)." );
  
  if ( _n_inserts_total == 0 ) {
    LogDate( "No inserts selected, so this assembly's HyperBasevector will be empty." );
    return true;
  }
  
  vecbasevector local_unibases;
  for ( int i = 0; i < _n_local_unipaths; i++ )
    local_unibases.push_back_reserve( _kbb.Seq( _hkp.EdgeObject( i ) ), 0, 2.0 );
  if (_DUMP_LOCAL_UNIBASES) {
    HyperBasevector hb( _hkp, _kbb );
    
    String hyper_file = _outhead + ".local_unibases.hbv";
    int fd = OpenForWrite( hyper_file );
    BinaryWrite( fd, hb );
    
    String unibases_file = _outhead + ".local_unibases.fasta";
    Ofstream( out, unibases_file );
    out << hb;

    String dot_file = _outhead + ".local_unibases.dot";
    Ofstream( outd, dot_file );
    _hkp.PrintSummaryDOT0w( outd, True, False, True, 0, True );
  }
  
  // Build an edited version of the HyperKmerPath from which hanging ends and
  // cyclic parts have been deleted.  This edited HKP will be merged into the
  // final assembly along with the insert walks.
  HyperKmerPath acyclic_HKP = MakeAcyclicHKP( _local_unipaths, local_unibases, min_size );
  LogDate( "Built direct acyclic assembly graph having " + 
	   ToString( acyclic_HKP.EdgeObjectCount( ) ) + " edges." );

  if ( _DUMP_LOCAL_UNIBASES ) {
    HyperBasevector hb( acyclic_HKP, _kbb );
    
    String hyper_file = _outhead + ".local_unibases_acyclic.hbv";
    int fd = OpenForWrite( hyper_file );
    BinaryWrite( fd, hb );
  }
  
  // Walk long inserts with fragments
  // The output HyperKmerPath is a union of all successfully walked inserts.
  // This is a timed operation!
  const int STRETCH = 3;
  double t_start = WallClockTime();
  vec<Bool> walked;
  HyperKmerPath inserts_HKP = WalkInserts( selected_pair_IDs, STRETCH, timer, _n_inserts_total, _n_inserts_walked, walked );

  _T_insert_walking = WallClockTime() - t_start;
  if (timer.TimedOut()) {
    // return false;
  }
  
  LogDate( ToString( _n_inserts_walked ) + " of " + ToString( _n_inserts_total )
	   + " (distinct) selected inserts successfully walked." );
  
  // Merge the walked inserts and the acyclic HKP into the HyperBasevector
  // object (hbv).
  // Note that the local HyperKmerPath (hkp) is NOT an input to MergeInsertWalksLG,
  // and its data will NOT necessarily be represented in hbv.  Only those portions
  // of hkp that span a walked insert (or come from the acyclic graph) will be
  // copied into hbv.  This keeps hbv free of extraneous, unverified data.
  // This is a timed operation!
  
  vec<HyperKmerPath> inserts_and_acyclic;
  if ( inserts_HKP.N() > 0 )
    inserts_and_acyclic.push_back( inserts_HKP );
  if (_USE_ACYCLIC && acyclic_HKP.N() > 0 ) inserts_and_acyclic.push_back(acyclic_HKP);
  t_start = WallClockTime();
  if ( inserts_and_acyclic.size() > 0 )
    MergeInsertWalks( inserts_and_acyclic, _local_unipaths, local_unibases );
  if (timer.TimedOut()) {
    // return false;
  }
  
  if ( _hbv.N() == 0 ) {
    LogDate( "HyperBasevector is empty - no local assembly here." );
    return true;
  }
  
  // If the merged HBV appears to contain duplication with individual components,
  // redo the computation without using acylic_HKP.  Horrid.
  
  if (_USE_ACYCLIC)
    {
    PerfectAlignerLG P( 1000, PerfectAlignerLG::findImproper );
    vec<alignment_plus> aligns;
    vecbasevector hbv_bases;
    for ( int i = 0; i < _hbv.EdgeObjectCount( ); i++ )
      hbv_bases.push_back_reserve( _hbv.EdgeObject(i) );
    P.Align( hbv_bases, aligns );
    if ( aligns.nonempty( ) ) { // perfect alignments -> duplication
      inserts_and_acyclic = vec<HyperKmerPath>( 1, inserts_HKP );
      MergeInsertWalks( inserts_and_acyclic, _local_unipaths, local_unibases );
    }
    }
  _T_insert_merging = WallClockTime() - t_start;
  
  _hbv.RemoveSmallComponents(min_size);
  
  // Call it done.
  
  LogDate( "Walked-insert neighborhoods merged into local assembly." );
  return true;
} 



/* Write the HyperBasevector of this assembly to files.
 *
 * INPUT: Filename
 * MUST ALREADY BE RUN: MakeInsertWalks
 *
 * OUTPUT: Writes the files <hbv_file> and <hbv_file>.fasta.  If LOCAL_DOT=True,
 * also writes <hbv_file>.dot.
 ******************************************************************************/
void
SeedNeighborhood::Write( size_t iSeedId, BinaryWriter& hbvWriter,
                         bool LOCAL_DOT, bool LOCAL_FASTA )
{
  LogDate( "Writing output files..." );
  
  // Write the HyperBasevector in binary.
  hbvWriter.write(static_cast<size_t>(iSeedId));
  hbvWriter.write(_hbv);

  // Write the HyperBasevector as a fastavector.
  if (LOCAL_FASTA)
  {    Ofstream( out, _outhead + ".hbv.fasta" );
       out << _hbv;    }
  
  if (LOCAL_DOT) {
    // Write the HyperBasevector as a DOT graph.
    Ofstream( out, _outhead + ".hbv.dot" );
    _hbv.PrintSummaryDOT0w( out, True, False, True, 0, True );
  }
  
  LogDate( "Finished logging output for seed ID " + ToString( _seed_ID ) );
  *_log << endl << endl;
}







/* Report on the unipath cloud - you must run MakeUnipathCloud before this!
 *
 * INPUT: Unipath placements on genome (optional - you can choose to load in
 *        an empty VecPlacementVec instead)
 *
 * REPORT:
 * 1. Seed unipath's length
 * 2. Number, total, size, and N50 size of the unipaths in this cloud
 * - if unipath_POGs was provided -
 * 3. Warning for unipaths that are not CN-1 or do not fall in this neighborhood
 * 4. Coverage of true neighborhood by unipath cloud
 * 5. List of unipaths, with predicted and actual location relative to the seed
 *
 ******************************************************************************/
void
SeedNeighborhood::EvalUnipathCloud( const vecbasevector & genome,
				    const VecPlacementVec& unipath_POGs ) const
{
  LogDate( "Running eval code: EvalUnipathCloud", 2 );
  ForceAssert( _n_cloud_unipaths != 0 );
  
  Bool evalReference = !unipath_POGs.empty() && !genome.empty();

  // Get the IDs of unipaths in this cloud
  vec<int> unipath_IDs( _n_cloud_unipaths );
  for ( int i = 0; i < _n_cloud_unipaths; i++ )
    unipath_IDs[i] = _cloud_unipaths[i].Uid( );
  
  // Tabulate lengths

  vec<int> lengths;
  lengths.reserve( _n_cloud_unipaths - 1 );
  for ( int i = 0; i < _n_cloud_unipaths; i++ ) {
    if ( unipath_IDs[i] == _seed_ID ) continue; // ignore the seed for now
    lengths.push_back( (*_global_unipath_lengths) [unipath_IDs[i]] );
  }
  
  
  // Report on length of seed and its cloud

  *_log << "\tSeed " << _seed_ID << " has length of "
       << (*_global_unipath_lengths) [_seed_ID] << " kmers." << endl;

  // The wording of the report changes if the cloud has 0, 1, or 2+ unipaths
  // aside from the seed

  switch ( lengths.isize( ) ) {
    
  case 0:
    *_log << "\tSeed's cloud contains no other unipaths." << endl;
    break;
    
  case 1:
    *_log << "\tSeed's cloud contains 1 other unipath, of size "
	 << lengths[0] << endl;
    break;
    
  default:
    Sort( lengths );
    *_log << "\tSeed's cloud contains " << _n_cloud_unipaths - 1
	 << " other unipaths; total length " << Sum( lengths )
	 << "; N50 length " << N50( lengths ) << endl;
  }
  
  int genome_ID = 0;
  int seed_loc = 0;
  vec<String> unipath_true_locs( _n_cloud_unipaths, "unknown" );

  // evaluations that are only possible with truth data

  if (evalReference) {
    
    *_log << "\n\tEvalUnipathCloud is examining unipath cloud's placement on reference..." << endl;
    
    // Report the seed unipath's actual copy number.
    PlacementVec::size_type seed_CN = unipath_POGs[_seed_ID].size( );
    *_log << "Seed " << _seed_ID << " has actual copy number " << seed_CN << "." << endl;
    
    // If the seed has any actual placements, use one of them to determine
    // the seed's "correct" location on the genome for reporting purposes.

    if ( seed_CN > 0 ) {
      genome_ID = unipath_POGs[_seed_ID][0].GenomeId( );
      seed_loc  = unipath_POGs[_seed_ID][0].pos( );
      
      // For each unipath in the cloud, write a tag describing its true location

      for ( int i = 0; i < _n_cloud_unipaths; i++ ) {
	PlacementVec places = unipath_POGs[ unipath_IDs[i] ];
	
	if ( places.empty( ) ) {
	  unipath_true_locs[i] = "Copy number 0";
	  continue;
	}
	
	// Find the placement of this unipath that's closest to the seed.

	placement closest_place = places[0];
	int place_dist = numeric_limits<int>::max( );
	for ( PlacementVec::size_type j = 0; j < places.size( ); j++ ) {
	  if ( places[j].GenomeId( ) != genome_ID ) continue;
	  int dist = Distance( places[j].Interval( ), unipath_POGs[_seed_ID][0].Interval( ) );
	  if ( dist < place_dist ) {
	    place_dist = dist;
	    closest_place = places[j];
	  }
	}
	
	unipath_true_locs[i] = ToString( closest_place.GenomeId( ) ) + "." 
	  + ToString( closest_place.pos( ) ) + "-" 
	  + ToString( closest_place.Pos( ) ) + " "
	  + ( closest_place.Fw( ) ? "fw" : "rc" ) + " " 
	  + ToString( places.size() );
      }
    }
  }
  
  
  // Write a chart of unipath locations (predicted and, if possible, true)
  
  *_log << endl;
  // Chart header line
  *_log << "\tUNIPATH seed neighborhood " << endl;
  *_log << "\tNote: the 'predicted loc' takes the true seed location as a reference." << endl << endl;
  *_log << "Local ID  Global ID   Size  predCN  Predicted loc                       True loc, orientation, trueCN" << endl;
  *_log << "--------  ---------  -----  ------  -------------                       -----------------------------" << endl;
  
  for ( int i = 0; i < _n_cloud_unipaths; i++ ) {
    int uID = _cloud_unipaths[i].Uid( );
    int start = _cloud_unipaths[i].Start( );
    int length = (*_global_unipath_lengths) [uID] + _K - 1;      
    int dev = _cloud_unipaths[i].MeanDev( );
    int CN = (*_global_predicted_CNs) [i];
    String seed_tag = ( uID == _seed_ID ) ? "SEED" : "    ";

    // Carefully format and print the chart row for this unipath    
    char predicted_loc[200];
    sprintf ( predicted_loc, "%d.%d-%d +/- %d",
	      genome_ID, start, ( start + length ), dev );
    char line[200];
    sprintf ( line, "%s %3d  %9d  %5d  %6d  %-35s %s", seed_tag.c_str( ), i, uID,
	      length, CN, predicted_loc, unipath_true_locs[i].c_str() );
    *_log << line << endl;
  }

  *_log << endl;
}



/* Report on the read cloud - you must run MakeReadCloud before this!
 *
 * INPUT: Unipath and read placements on genome (optional - you can choose to
 *        load in an empty VecPlacementVec and vec<ReadLocationLG> instead)
 *
 * REPORT:
 * 1. Number of reads in cloud
 * 2. Coverage of seed unipath by the read cloud
 * - if read_POGs and unipath_POGs was provided -
 * 3. Number of reads that fall outside of this neighborhood
 * 4. Coverage of true neighborhood by seed cloud
 *
 ******************************************************************************/
void
SeedNeighborhood::EvalReadCloud( const vecbasevector & genome, 
				 const VecPlacementVec& unipath_POGs,
				 const vec<ReadLocationLG> & read_POGs ) const
{
  LogDate( "Running eval code: EvalReadCloud", 2 );
  ForceAssert( _n_cloud_unipaths != 0 );
  if ( _n_reads == 0 ) {
    *_log << "\tThe local read cloud is empty.  No evaluation to be done.  (Run MakeReadCloud?)" << endl;
    return;
  }
  
  *_log << "\t" << _n_reads << " reads in cloud." << endl;
  
  // Make an index of local read ID -> ReadLocationLG ID in unilocs vector
  vec<int> read_to_uniloc( _n_reads, -1 );
  for ( int i = 0; i < (*_global_unilocs).isize( ); i++ ) {
    longlong pos = BinPosition( _local_to_global, (*_global_unilocs)[i].ReadId( ) );
    if ( pos == -1 ) continue; // i.e. this read isn't in this nhoodx
    read_to_uniloc[pos] = i;
  }
  
  // Examine the unilocs that contain a read in this neighborhood
  vec<longlong> reads_on_seed; // list of reads that cover the seed unipath
  reads_on_seed.reserve( _n_reads );
  for ( int i = 0; i < _n_reads; i++ ) {
    int uniloc_id = read_to_uniloc[i];
    ReadLocationLG uniloc = (*_global_unilocs)[ uniloc_id ];
    if ( uniloc.Contig( ) == _seed_ID )
      reads_on_seed.push_back( i );
  }
  
  
  // Examine coverage of the seed unipath (note: in base space, not kmer space)
  int seed_length = (*_global_unipath_lengths)[_seed_ID] + _K - 1;
  vec<unsigned short> seed_cov( seed_length, 0 );
  
  for ( int i = 0; i < reads_on_seed.isize( ); i++ ) {
    ReadLocationLG uniloc = (*_global_unilocs)[ read_to_uniloc[ reads_on_seed[i] ] ];
    if ( uniloc.Contig( ) != _seed_ID ) continue; // this should never happen
    
    // Mark this portion of the seed unipath as covered
    int start = uniloc.StartOnContig( );
    int stop = start + (*_global_read_lengths)[_local_to_global[reads_on_seed[i]]];
    if ( start < 0 ) start = 0;
    if ( stop > seed_length ) stop = seed_length;
    
    for ( int j = start; j < stop; j++ )
      if ( seed_cov[j] + 1 != 0 ) // avoid overflow
	seed_cov[j]++;
  }
  
  int n_holes = seed_cov.CountValue( 0 );
  
  // Report on seed unipath's coverage
  *_log << "\tSeed unipath " << _seed_ID << " aligns to "
       << reads_on_seed.isize( ) << " reads in its cloud." << endl;
  *_log << "\tThese reads give the seed an average base coverage of "
       << ( (double)Sum( seed_cov ) ) / seed_length << "." << endl;
  if ( n_holes > 0 )
    *_log << "\t" << n_holes << " of " << seed_length
	 << " bases on the seed are not covered by the read cloud!" << endl;
  else
    *_log << "\tAll " << seed_length
	 << " bases on the seed are covered by the read cloud." << endl;
  
  
  
  
  // The rest of this analysis requires read and unipath placements on genome
  if ( !      genome.size( ) ) return;
  if ( !unipath_POGs.size( ) ) return;
  if ( !   read_POGs.size( ) ) return;
  
  *_log << "\n\tEvalReadCloud is examining read cloud's placement on reference..." << endl;
  
  // Report the seed unipath's actual copy number.
  PlacementVec::size_type seed_CN = unipath_POGs[_seed_ID].size( );
  *_log << "Seed " << _seed_ID << " has actual copy number " << seed_CN << "." << endl;
  if ( seed_CN != 1 ) return; // The remaining analysis won't work for a non-CN1 seed
  
  placement const& seed_place = unipath_POGs[_seed_ID][0];
  int genome_ID = seed_place.GenomeId( );
  
  // Re-order the read locations to suit the read IDs in this cloud
  vec<vec<const ReadLocationLG *> > local_read_POGs( _n_reads );
  for ( int i = 0; i < read_POGs.isize( ); i++ ) {
    longlong pos = BinPosition( _local_to_global, read_POGs[i].ReadId( ) );
    if ( pos == -1 ) continue; // i.e. this read isn't in this nhood
    local_read_POGs[pos].push_back( &read_POGs[i] );
  }
  
  
  // Set up vector to mark coverage of this neighborhood by its reads
  int nhood_offset = Max( seed_place.pos( ) - _NHOOD_RADIUS_INTERNAL, 0 );
  int nhood_max    = Min( seed_place.Pos( ) + _NHOOD_RADIUS_INTERNAL,
			  genome[genome_ID].isize( ) );
  int nhood_size = nhood_max - nhood_offset;
  vec<bool> nhood_cov( nhood_size, false );
  
  // Appraise each read in this cloud, and calculate neighborhood coverage
  int n_reads_high_CN = 0;
  int n_reads_on_wrong_contig = 0;
  int n_reads_outside_radius = 0;
  int n_good_reads = 0;
  
  
  for ( int i = 0; i < _n_reads; i++) {
    const vec<const ReadLocationLG *> locs = local_read_POGs[i];
    
    
    // Ignore high-CN reads (they might be fine)
    if ( locs.isize( ) != 1 ) {
      n_reads_high_CN++;
      continue;
    }
    
    // Check that the read is on the same contig of genome as the seed
    if ( locs[0]->Contig( ) != genome_ID ) {
      n_reads_on_wrong_contig++;
      continue;
    }
    
    // Check that the read's placement is within the neighborhood radius
    int read_length = (*_global_read_lengths)[ _local_to_global[i] ];
    int start = locs[0]->StartOnContig( );
    int dist_from_seed;
    if ( start > seed_place.Pos( ) )
      dist_from_seed = start - seed_place.Pos( ) - _K + 1;
    else
      dist_from_seed = seed_place.pos( ) - start - read_length;
    
    static const int link_size = 4000; // TODO: de-hardwire this
    if ( dist_from_seed > _NHOOD_RADIUS_INTERNAL + link_size ) {
      n_reads_outside_radius++;
      continue;
    }
    
    // This is a good read!  Mark its region of the neighborhood as covered
    int start_in_nhood = start - nhood_offset;
    int stop_in_nhood  = start_in_nhood + read_length;
    if ( start_in_nhood < 0 ) start_in_nhood = 0;
    if ( stop_in_nhood > nhood_size ) stop_in_nhood = nhood_size;
    for ( int j = start_in_nhood; j < stop_in_nhood; j++ )
      nhood_cov[j] = true;
    
    n_good_reads++;
  }
  
  // Report on good/bad reads in this cloud
  double pct_good_reads = 100.0 * ( (double) n_good_reads / _n_reads );
  *_log << "\tTally of bad reads in this cloud:" << endl;
  *_log << "\t(" << n_reads_high_CN
       << " high-copy-number reads ignored for this analysis.)" << endl;
  *_log << "\tReads that belong on a different reference contig: "
       << n_reads_on_wrong_contig << endl;
  *_log << "\tReads that belong outside the radius (NHOOD_RADIUS_INTERNAL = "
       << _NHOOD_RADIUS_INTERNAL << ") of the seed unipath: "
       << n_reads_outside_radius << endl;
  *_log << "\tThe remaining " << n_good_reads << " reads (" << pct_good_reads
       << "%) are in the neighborhood." << endl;
  
  // Report on neighborhood coverage by the read cloud
  double nhood_coverage = 100.0 * ( (double) Sum( nhood_cov ) / nhood_size );
  *_log << "\n\tNeighborhood size, for reads, is " << nhood_size
       << " (measured as seed size of " << seed_place.Pos( ) - seed_place.pos( )
       << " bases,\n\t\tplus an extension of up to NHOOD_RADIUS_INTERNAL = "
       << _NHOOD_RADIUS << " on either side)" << endl;
  *_log << "\tNeighborhood bases are " << nhood_coverage
       << "% covered by the read cloud." << endl;
}


/* Report on the local assembly - you must run MakeLocalAssembly before this!
 *
 * REPORT:
 * 1. Lengths of unipaths in local assembly HyperKmerPath
 * 2. Sizes of connected components in local assembly HyperKmerPath
 *
 ******************************************************************************/
void
SeedNeighborhood::EvalLocalAssembly( ) const
{
  LogDate( "Running eval code: EvalLocalAssembly", 2 );
  if ( _n_local_unipaths == 0 ) {
    *_log << "\tThere are no local unipaths.  No evaluation to be done.  (Run MakeLocalAssembly?)" << endl;
    return;
  }
  
  // Find lengths of local unipaths
  vec<int> local_lengths( _n_local_unipaths, 0 );
  for ( int i = 0; i < _n_local_unipaths; i++ )
    local_lengths[i] = _hkp.EdgeLength( i );
  
  // Report on unipath lengths
  *_log << "\t\t\t\tLocal assembly HyperKmerPath contains " << _n_local_unipaths
       << " unipaths; longest is " << Max( local_lengths )
       << " kmers " << endl;
  
  *_log << "Unipath lengths (kmers):" << endl;
  PrintBasicNStats( "", local_lengths, *_log );
  
  
  // Report unipath graph connected component sizes
  vec< vec<int> > comps;
  _AG.Components( comps );  // get components
  vec<int> compSizes( comps.isize(), 0 ); // component sizes
  int unipathsSize = 0; // sum unipath sizes
  for ( int cid = 0; cid < comps.isize(); cid++ ){
    long compSize = 0;
    for ( int j = 0; j < comps[ cid ].isize(); j++ )
      compSize += local_lengths[ comps[ cid ][ j ] ];

    unipathsSize      += compSize;
    compSizes[ cid ]  = compSize;
  }
  ReverseSort( compSizes );
  
  *_log << "\nUnipath adjacency graph contains " << comps.isize() << " components " << endl;
  *_log << "Largest component sizes (those with more than 1.0% of all component kmers):" << endl;
  for ( int i = 0; i < compSizes.isize(); i++ ){
    double percSize = 100.0 * (double) compSizes[ i ] / (double) unipathsSize;
    if( percSize < 1.0 )  continue;
    *_log << "\t" << comps[ i ].isize() << " unipaths,  " << compSizes[i] << " kmers,  " << percSize << " % of all kmers" << endl;
  }
  
  // Report on component sizes
  *_log << "\nHyperKmerPath component sizes (i.e., kmers per connected component):" << endl;
  PrintBasicNStats( "", compSizes, *_log );

  // Find lengths of cloud unipaths (for comparison)
  vec<int> cloud_lengths( _n_cloud_unipaths, 0 );
  for ( int i = 0; i < _n_cloud_unipaths; i++ )
    cloud_lengths[i] = (*_global_unipath_lengths) [ _cloud_unipaths[i].Uid( ) ];
  
  Sort( compSizes );
  *_log << "\t\t\t\tN50 component length in local assembly: "
       << N50( compSizes ) << endl;
  *_log << "\t\t\t\t  ... as compared to in cloud unipaths: "
       << N50( cloud_lengths ) << endl;
}



/* Report on the local insert-walks - you must run MakeInsertWalks before this!
 *
 * REPORT:
 * 1. Lengths of unipaths in local assembly HyperBasevector
 * 2. Sizes of connected components in local assembly HyperBasevector
 *
 ******************************************************************************/
void
SeedNeighborhood::EvalInsertWalks( ) const
{
  LogDate( "Running eval code: EvalInsertWalks", 2 );
  if ( _hbv.N( ) == 0 ) {
    *_log << "\tThere is no local HyperBasevector.  No evaluation to be done.  (Run MakeInsertWalks?)" << endl;
    return;
  }
  
  // Find lengths of bases in local HyperBasevector
  vec<basevector> bases = _hbv.Edges( );
  int n_edges = bases.isize( );
  vec<int> base_lengths;
  base_lengths.reserve( n_edges );
  for ( int i = 0; i < n_edges; i++ )
    base_lengths.push_back( bases[i].size( ) );
  
  // Report on unipath lengths
  *_log << "\t\t\t\tLocal assembly HyperBasevector contains " << n_edges
       << " unibases; longest is " << Max( base_lengths ) << " bases " << endl;
  
  *_log << "Basevector lengths (bases):" << endl;
  PrintBasicNStats( "", base_lengths, *_log );
  
  
  // Find component sizes in HyperBasevector
  // Note the two ways to measure the "size" of a component: the number of edges
  // (basevectors) in it, and the total length of all the bases in it.
  vec<vec<int> > components;
  _hbv.ComponentEdges( components );
  int n_components = components.isize( );
  vec<int> component_sizes_edges( n_components );
  vec<int> component_sizes_bases( n_components );
  for ( int i = 0; i < n_components; i++ ) {
    component_sizes_edges[i] = components[i].isize( );
    int size = 0;
    for ( int j = 0; j < components[i].isize( ); j++ )
      size += _hbv.EdgeObject( components[i][j] ).size( );
    component_sizes_bases[i] = size;
  }
  
  // Report on component sizes (both forms)
  *_log << "\nHyperBasevector component sizes (EDGES per connected component):" << endl;
  PrintBasicNStats( "", component_sizes_edges, *_log );
  *_log << "\nHyperBasevector component sizes (BASES per connected component):" << endl;
  PrintBasicNStats( "", component_sizes_bases, *_log );
}





// Helper function: Print to the logstream with a date tag
void
SeedNeighborhood::LogDate( const String & s, const int n_newlines ) const
{
  for ( int i = 0; i < n_newlines; i++ ) *_log << endl;
  *_log << Date( ) << ":\t" << s << endl;
}
