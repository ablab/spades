///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// RecoverUnipaths.  Recover some of the unipaths that do not occur in the
// final assembly.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "math/Functions.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"



// MakeDepend: dependency QueryLookupTable

Bool cmp_len_rev( const basevector& b1, const basevector& b2 )
{    if ( b1.size( ) > b2.size( ) ) return True;
     if ( b1.size( ) < b2.size( ) ) return False;
     return b1 > b2;    }



vec<int>
RemoveUsedUnipaths( const String & reads_head, const String & sub_dir, const int K, const vecKmerPath & paths, const vecKmerPath & paths_rc, const vec<tagged_rpint> & pathsdb )
{
  vec<int> used_to_delete;
  vec<int> to_rc;
  size_t unibase_count; 
  vecbitvector cov;
  
  {
    vecbasevector unibases( reads_head + ".unibases.k" + ToString(K) );
    unibase_count = unibases.size();
    UnibaseInvolution( unibases, to_rc, K );
    Mimic( unibases, cov );
  }
  
  cout << Date( ) << ": unibases has size " << unibase_count << endl;
  
  size_t n_hyper_paths = MastervecFileObjectCount( reads_head + ".paths.k" + ToString(K) );
  size_t orig_hyper_paths = n_hyper_paths - unibase_count;
  
  HyperKmerPath h( sub_dir + "/hyper" );
  
  // Look for assembly kmers in the unipaths
  for ( int i = 0; i < h.EdgeObjectCount( ); i++ ) {
    const KmerPath& p = h.EdgeObject(i);
    
    for ( int j = 0; j < p.NSegments( ); j++ ) {
      vec<longlong> con;
      Contains( pathsdb, p.Segment(j), con );
      
      for ( int l = 0; l < con.isize( ); l++ ) {
	const tagged_rpint& t = pathsdb[ con[l] ];
	int path_id = t.ReadId( );
	int unipath_id = t.ReadId( ) - orig_hyper_paths;
	
	// Is this a unipath
	if ( unipath_id < 0 ) continue;
	
	// Obtain the unipath
	const KmerPath & kp = (t.Fw() ? paths[path_id] : paths_rc[path_id]);
	
	// Find the starting kmer in the unipath of this segment
	size_t interval_start = 0;
	for (size_t interval = 0; interval < t.PathPos(); interval++) {
	  interval_start += kp.Length(interval);
	}
	
	// Add coverage
	for ( longlong r = p.Segment(j).Start( ); r <= p.Segment(j).Stop( ); r++ ) {
	  if ( r >= t.Start( ) && r <= t.Stop( ) ) {
	    cov[unipath_id].Set( r - t.Start( ) + interval_start, True );
	  }
	} 
      } 
    } 
  }
  
  for ( size_t i = 0; i < unibase_count; i++ )
    if ( Sum(cov[i]) > ( cov[i].size( ) - K + 1 ) / 2 )
      used_to_delete.push_back( i, to_rc[i] );
  
  UniqueSort(used_to_delete);
  cout << Date( ) << ": used_to_delete has size " << used_to_delete.size( )
       << endl;
  return used_to_delete;
}






void LongPaths( const HyperKmerPath& h, const int min_len, 
     vec<KmerPath>& LONG_PATHS )
{    
     // Set up initial list to process.

     vec<HyperKmerPath> to_process;
     to_process.push_back(h);

     // Process each HyperKmerPath.

     while( to_process.nonempty( ) )
     {    HyperKmerPath hi = to_process.back( );
          to_process.pop_back( );

          // Find long paths.

          vec<int> D;
          const int infinity = 2000000000;
          DistancesToEnd( hi, &KmerPath::KmerCount, infinity, True, D );
          vec< pair<int,int> > dv;
          for ( int v = 0; v < hi.N( ); v++ )
               dv.push( D[v], v );
          ReverseSort(dv);
          if ( dv[0].first < min_len ) continue;
          vec< vec<int> > long_paths;
          set<int> used;
          for ( size_t i = 0; i < dv.size(); i++ )
          {    if ( dv[i].first < min_len ) break;
               int v = dv[i].second;
               vec<int> long_path;
               Bool fail = False;
               while( hi.From(v).nonempty( ) )
               {    for ( size_t j = 0; j < hi.From(v).size(); j++ )
                    {    int w = hi.From(v)[j];
                         int e = hi.EdgeObjectIndexByIndexFrom(v, j);
                         if ( D[w] == D[v] - hi.EdgeObject(e).KmerCount( ) )
                         {    if ( Member( used, e ) )
                              {    fail = True;
                                   break;    }
                              long_path.push_back(e);
                              v = w;
                              break;    }    }
                    if (fail) break;    }
               if ( !fail )
               {    long_paths.push_back(long_path);
                    for ( size_t j = 0; j < long_path.size(); j++ )
                         used.insert( long_path[j] );    }    }

          // Make them into KmerPaths.

          vec<KmerPath> Ls;
          for ( size_t i = 0; i < long_paths.size(); i++ )
          {    KmerPath L;
               for ( size_t j = 0; j < long_paths[i].size(); j++ )
                    L.Append( hi.EdgeObject( long_paths[i][j] ) );
               Ls.push_back(L);    }
          LONG_PATHS.append(Ls);

          // Find edges that are part of an alternate path from the beginning
          // of long path to its end.  These will be treated as part of
          // the paths and deleted.

          vec<int> to_left, to_right, suc, pre;
          hi.ToLeft(to_left), hi.ToRight(to_right);
          vec<int> stuff;
          for ( size_t j = 0; j < long_paths.size(); j++ )
          {    hi.GetSuccessors1( to_left[ long_paths[j].front( ) ], suc );
               hi.GetPredecessors1( to_right[ long_paths[j].back( ) ], pre );
               vec<int> between = Intersection( suc, pre );
               for ( size_t i = 0; i < between.size(); i++ )
               {    int v = between[i];
                    for ( size_t r = 0; r < hi.From(v).size(); r++ )
                    {    int w = hi.From(v)[r];
                         if ( BinMember( between, w ) )
                         {    stuff.push_back(
                                   hi.EdgeObjectIndexByIndexFrom( v, r ) );
                                        }    }    }    }
          UniqueSort(stuff);
               
          // Delete edges in long paths and their reverse complements.

          vec<KmerPath> edges;
          for ( size_t l = 0; l < stuff.size(); l++ )
               edges.push_back( hi.EdgeObject( stuff[l] ) );
          UniqueSort(edges);
          vec<int> to_delete(stuff);
          for ( int l = 0; l < hi.EdgeObjectCount( ); l++ )
          {    if ( BinMember( stuff, l ) ) continue;
               KmerPath p = hi.EdgeObject(l);
               p.Reverse( );
               if ( BinMember( edges, p ) ) to_delete.push_back(l);    }
          hi.DeleteEdges(to_delete);
          hi.RemoveDeadEdgeObjects( );

          // Find components in what's left and save.

          vec< vec<int> > componentsi;
          hi.Components(componentsi);
          for ( size_t j = 0; j < componentsi.size(); j++ )
          {    HyperKmerPath hij( hi, componentsi[j] );
               to_process.push_back(hij);    }    }    }



vecbasevector
FindRecoverableSeqs( const HyperKmerPath & hkp, const KmerBaseBroker & kbb, const int MIN_KEEPER )
{
  cout << Date( ) << ": linearizing graph" << endl;
  vec< vec<int> > components;
  hkp.Components(components);
  size_t n_components = components.size();
  const int slen = 40;
  cout << Date( ) << ": " << n_components << " components" << endl;
  // Bases/paths in each component
  vec<vecbasevector> c_bases( n_components );
  
#pragma omp parallel for
  for ( size_t ci = 0; ci < n_components; ci++ )
    {    
      // Find long paths.
      
      HyperKmerPath hi( hkp, components[ci] );
      vec<basevector> long_seqs0;
      {
	vec<KmerPath> long_paths0;
	LongPaths( hi, MIN_KEEPER - hkp.K() + 1, long_paths0 );
	for ( size_t i = 0; i < long_paths0.size(); i++ )
	  long_seqs0.push_back( kbb.Seq( long_paths0[i] ) );
      }
      sort( long_seqs0.begin(),  long_seqs0.end(),  cmp_len_rev );
      
      // Identify redundancy.
      
      vec< triple< basevector, size_t, size_t > > X;
      for ( size_t i = 0; i < long_seqs0.size(); i++ )
	{    for ( size_t j = 0; j <= long_seqs0[i].size() - slen; j++ )
	    {    basevector b;
	      b.SetToSubOf( long_seqs0[i], j, slen );
	      X.push( b, i, j );
	      b.ReverseComplement( );
	      X.push( b, i, j );    }    }
      Sort(X);
      vec< vec<ho_interval> > dels( long_seqs0.size( ) );
      for ( size_t i = 0; i < X.size(); i++ )
	{    size_t j;
	  for ( j = i + 1; j < X.size(); j++ )
	    if ( X[j].first != X[i].first ) break;
	  for ( size_t k = i + 1; k < j; k++ )
	    {    if ( X[k].second > X[i].second )
		{    dels[ X[k].second ].push( 
					      X[k].third, X[k].third + slen );    }    }
	  i = j - 1;    }
      
      // Delete dinucleotide 40-mers.
      
      vec<basevector> dinukes;
      for ( int i = 0; i < 4; i++ )
	{    for ( int j = 0; j < 4; j++ )
	    {    if ( j == i ) continue;
	      basevector b(slen);
	      for ( int r = 0; r < slen; r++ )
		{    if ( r % 2 == 0 ) b.Set( r, i );
		  else b.Set( r, j );    }
	      dinukes.push_back(b);    }    }
      Sort(dinukes);
      for ( size_t i = 0; i < X.size(); i++ )
	{    if ( BinMember( dinukes, X[i].first ) )
	    {    dels[ X[i].second ].push( X[i].third, X[i].third + slen );    }    }
      
      // Delete bad stuff.
      
      for ( size_t i = 0; i < long_seqs0.size(); i++ )
	{
	  vec<ho_interval> un;
	  Uncovered( long_seqs0[i].size( ), dels[i], un );
	  for ( size_t j = 0; j < un.size(); j++ )
	    {    if ( un[j].Length( ) >= MIN_KEEPER )
		{
		  // Find the not-bad basevector, and the corresponding not-bad KmerPath.
		  basevector b;
		  b.SetToSubOf( long_seqs0[i], un[j].Start( ), un[j].Length( ) );
		  c_bases[ci].push_back_reserve(b);
		}    }    }    }
  
  // Combine all the linear sequences.
  vecbasevector seqs;
  for ( size_t ci = 0; ci < n_components; ci++ )
    seqs.Append( c_bases[ci] );
  return seqs;
}




int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String(SUBDIR);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(WRUN_IN, "");
     CommandArgument_String_OrDefault(WRUN_OUT, "recover");
     CommandArgument_String_OrDefault( READS, "reads" );
     CommandArgument_String_OrDefault( HYPER_OUT, "hyper_plus" );
     CommandArgument_Bool_OrDefault(USE_TRUTH, False);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Int_OrDefault(MIN_KEEPER, 1000);
     CommandArgument_Bool_OrDefault( MEM_DEBUG, False );
     EndCommandArguments;

     // Thread control
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String KS = ToString(K);
     String reads_head = sub_dir + "/" + READS;
     String refname = data_dir + "/genome";

     // Load assembly, etc.

     if ( MEM_DEBUG ) PrintDateAndMemUsage();
     
     cout << Date() << ": Loading reads paths" << endl;
     vecKmerPath paths   ( reads_head + ".paths.k" + KS );
     vecKmerPath paths_rc( reads_head + ".paths_rc.k" + KS );
     BREAD2( reads_head + ".pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
     ForceAssertEq( paths.size(), paths_rc.size() );
          
     // Build unipath adjacency hyperkmerpath graph.
     // This is the object we will be doing most of our work on:
     // we will remove unipaths from it as they fail to meet the criteria for
     // recoverability, and then recover the few remaining ones.

     vecbasevector recovered_seqs;
     {    HyperKmerPath hkp;
          {    vecKmerPath unipaths( reads_head + ".unipaths.k" + ToString(K) );
               digraph A;
               BinaryRead( reads_head + ".unipath_adjgraph.k" + ToString(K), A );
               cout << Date( ) << ": building unipath adjacency HyperKmerPath" 
                    << endl;
               BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, hkp );    }
          if ( MEM_DEBUG ) PrintDateAndMemUsage();

          // Remove used edges.  An edge of length n kmers is declared used if the 
          // number of its kmers that appear in the assembly exceeds n/2.

          cout << Date( ) << ": Removing used unipaths" << endl; 
          {    vec<int> used_to_delete = RemoveUsedUnipaths( reads_head, sub_dir, K, 
                    paths, paths_rc, pathsdb );
               hkp.DeleteEdges(used_to_delete);    }
          hkp.RemoveDeadEdgeObjects( );
          if ( MEM_DEBUG ) PrintDateAndMemUsage();

          // Remove the loop subgraph (e.g., the subgraph of this HKP consisting of
          // all internal cycles.)  We forcibly add the reverse complement of
          // loops, which shouldn't be necessary and suggests that there is a bug in
          // LoopEdges.

          cout << Date( ) << ": Removing loop subgraph" << endl;
          hkp.RemoveLoopSubgraph();

          // Clean up.
     
          hkp.RemoveUnneededVertices( );
          hkp.RemoveDeadEdgeObjects( );
          hkp.RemoveEdgelessVertices( );
          if ( MEM_DEBUG ) PrintDateAndMemUsage();

          // Delete any palindromic edges.

          cout << Date( ) << ": Removing palindromic edges" << endl;
          {    vec<int> palindromic;
               for ( int i = 0; i < hkp.EdgeObjectCount( ); i++ )
               {    KmerPath p = hkp.EdgeObject(i), q;
                    p.Canonicalize( );
                    q = p;
                    q.Reverse( );
                    q.Canonicalize( );
                    if ( p == q ) palindromic.push_back(i);    }
               hkp.DeleteEdges(palindromic);    }
          hkp.RemoveDeadEdgeObjects( );
          hkp.RemoveEdgelessVertices( );
          if ( MEM_DEBUG ) PrintDateAndMemUsage();

          // If two components of the graph are reverse complements of each other,
          // delete one of them.

          hkp.DeleteReverseComplementComponents( );

          // Remove hanging ends and small components.  We employ a particularly 
          // stringent test for component size.

          cout << Date( ) << ": Removing hanging ends" << endl;
          RemoveHangingEnds( hkp, &KmerPath::KmerCount, 250, 5.0 );
          const int min_component = MIN_KEEPER - K + 1;
          hkp.RemoveSmallComponents(min_component);
          {    vec<int> dist_to_end;
	       DistancesToEnd( hkp, &KmerPath::KmerCount, min_component, True, 
                    dist_to_end );
               vec< vec<int> > components;
               hkp.Components(components);
               vec<int> edges_to_delete;
               for ( size_t i = 0; i < components.size(); i++ )
               {    int D = 0;
                    for ( size_t j = 0; j < components[i].size(); j++ )
                         D = Max( D, dist_to_end[ components[i][j] ] );
                    if ( D < min_component )
                    {    edges_to_delete.append( 
                              hkp.EdgesBetween( components[i] ) );    }    }
               hkp.DeleteEdges(edges_to_delete);     }
          hkp.RemoveUnneededVertices( );
          hkp.RemoveDeadEdgeObjects( );
          if ( MEM_DEBUG ) PrintDateAndMemUsage();
     
          // Find the list of unipaths to recover.

          cout << Date( ) << ": Generating a list of recoverable unipaths" << endl;
          {    KmerBaseBroker kbb( K, paths, paths_rc, pathsdb, 
                    reads_head + ".fastb" );
               recovered_seqs = FindRecoverableSeqs( 
                    hkp, kbb, MIN_KEEPER );    }    }
     
     // Destroy some no-longer-needed structures.

     Destroy( paths );
     Destroy( paths_rc );
     Destroy( pathsdb );
     if ( MEM_DEBUG ) PrintDateAndMemUsage();
     
     // Write the sequences.

     size_t n_recovered_seqs = recovered_seqs.size();
     cout << Date( ) << ": found " << n_recovered_seqs << " sequences of total size "
          << recovered_seqs.sumSizes() << endl;
     String outname = ( WRITE ? "hyper.extra.fastb" : "hyper.extra.test.fastb" );
     recovered_seqs.WriteAll( sub_dir + "/" + outname );
     if ( MEM_DEBUG ) PrintDateAndMemUsage();

     // Build the output HyperKmerPath.
     // The output HyperKmerPath is a combination of the assembly made by
     // MergeNeighborhoods and the sequences we've recovered from the unipaths.

     if (WRITE)
     {    String wrun_in_dir = sub_dir + "/" + WRUN_IN;
          vec<HyperBasevector> h2(2);
          {    vec<basevector> seqsx;
               for ( size_t i = 0; i < recovered_seqs.size( ); i++ )
                    seqsx.push_back( recovered_seqs[i] );
               KmerBaseBroker kbb( wrun_in_dir, K );
               HyperKmerPath h( sub_dir + "/hyper" );
               ForceAssertEq( K, h.K( ) );
               h2[0] = HyperBasevector( h, kbb ); 
               h2[1] = HyperBasevector( K, seqsx );    }
          HyperBasevector hbplus( K, h2 );
          Destroy(h2);
	  String wrun_out_dir = sub_dir + "/" + WRUN_OUT;
          Mkdir777(wrun_out_dir);
	  // write to wrun_out_dir so we can use them again in kbb2.
          vecKmerPath spaths;
	  if ( MEM_DEBUG ) PrintDateAndMemUsage();
          {    vecbasevector bases;
               bases.reserve( bases.size( ) + hbplus.EdgeObjectCount( ) );
               for ( int i = 0; i < hbplus.EdgeObjectCount( ); i++ )
                    bases.push_back_reserve( hbplus.EdgeObject(i) );
               bases.WriteAll( wrun_out_dir + "/reads.fastb" );
	       // Re-path (runtime bottleneck)
	       cout << Date() << ": Re-pathing (for later use in KBBs)" << endl;
               ReadsToPathsCoreY( bases, K, spaths, sub_dir + "/tmp",
                    omp_get_max_threads( ) );    }
	  if ( MEM_DEBUG ) PrintDateAndMemUsage();
          spaths.WriteAll( wrun_out_dir + "/reads.paths.k" + KS );
          {    vecKmerPath spaths_rc(spaths);
               for ( size_t i = 0; i < spaths.size( ); i++ )
                    spaths_rc[i].Reverse( );
               spaths_rc.WriteAll( wrun_out_dir + "/reads.paths_rc.k" + KS );
               vec<tagged_rpint> spathsdb;
               CreateDatabase( spaths, spaths_rc, spathsdb );
               BinaryWrite2( wrun_out_dir + "/reads.pathsdb.k" + KS, spathsdb );    }
	  
	  cout << Date( ) << ": writing " << HYPER_OUT << " and auxiliary files" 
               << endl;
	  vec<KmerPath> pathsx = VecOfKmerPath( spaths );
	  HyperKmerPath hplus( K, hbplus, pathsx );
	  String out_head = sub_dir + "/" + HYPER_OUT;
	  BinaryOverwrite( out_head, hplus );
	  
	  KmerBaseBroker kbb2( wrun_out_dir, K );
	  hplus.DumpFasta( out_head + ".fasta", kbb2 );
	  hplus.DumpFastb( out_head + ".fastb", kbb2 );
	  hplus.DumpGraphML( out_head + ".graphml" );
     }
     if ( MEM_DEBUG ) PrintDateAndMemUsage();

     // Align to reference.
     
     if ( USE_TRUTH && IsRegularFile( refname + ".lookup" ) ) {
       if (recovered_seqs.empty( ))
	 System( "touch " + sub_dir + "/hyper.extra.to_align" ); // keep MakeMgr happy
       else {
	 cout << Date( ) << ": aligning to reference (optional)" << endl;
	 vec<look_align> aligns;
	 PerfectLookup( 12, recovered_seqs, refname + ".lookup", aligns, FW_OR_RC );
	 vec<Bool> aligned( n_recovered_seqs, False );
	 for ( size_t i = 0; i < aligns.size(); i++ )
	   aligned[ aligns[i].query_id ] = True;
	 Ofstream( tout, sub_dir + "/hyper.extra.to_align" );
	 for ( size_t i = 0; i < n_recovered_seqs; i++ )
	     if ( !aligned[i] ) tout << i << "\n"; 
	 cout << "\nALIGNMENT OF IMPERFECT SEQUENCES TO REFERENCE:\n";
	 SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.05 SEQS=" + sub_dir +
			+ "/" + outname + " L=" + refname + ".lookup VISUAL=True "
			+ "LIST_UNPLACED=True " + "TMP_DIR=" + sub_dir + "/tmp "
			+ "ALIGN_UNALIGNED_BITS=True MF=500:500 QUIET=True NH=True "
			+ "SEQS_TO_PROCESS=@" + sub_dir+ "/hyper.extra.to_align" );
       } 
     }
     if ( MEM_DEBUG ) PrintDateAndMemUsage();
     
     cout << Date( ) << ": done!" << endl;   
}
