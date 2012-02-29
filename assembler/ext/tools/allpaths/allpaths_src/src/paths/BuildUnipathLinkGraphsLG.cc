///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// This code builds the unipath link graphs used by LocalizeReads.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "MainTools.h"
#include "ReadLocationLG.h"
#include "Set.h"
#include "graph/FindCells.h"
#include "paths/GetNexts.h"
#include "paths/Unipath.h"
#include "paths/UnipathNhoodLG.h"

template<int K0>
void GetKmers( const vecbasevector& unibases, vec< kmer_record<K0,2> >& kmers )
{    kmers.clear( );
     for ( ulonglong i = 0; i < unibases.size( ); i++ )
     {    for ( ulonglong j = 0; j <= unibases[i].size( ) - K0; j++ )
          {    basevector b;
               b.SetToSubOf( unibases[i], j, K0 );
               kmer_record<K0,2> rec;
               rec.Set( b, i, j );
               kmers.push_back(rec);    }    }
     Sort(kmers);    }

int main( int argc, char** argv ) {
  
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_String_OrDefault_Doc( GRAPH, "unipath_link_graph", "Prefix for file name in SUBDIR" );
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault(MIN_EDGE_MULTIPLICITY, 3);
  CommandArgument_Bool_OrDefault(TRANSITIVE_FILL_IN, True);
  CommandArgument_Int_OrDefault(TRANSITIVE_RADIUS, 6000);
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Bool_OrDefault(PRINT_LINKS, False);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Int_OrDefault(MIN_UNIPATH_FOR_TRANSITIVE_JOIN, 0);
  CommandArgument_Bool_OrDefault(USE_ADJACENCIES_AS_LINKS, True);
  CommandArgument_Int_OrDefault(NORMAL_MULT, 2);
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = PRE + "/" + DATA + "/" + RUN + "/ASSEMBLIES/" + SUBDIR;
  String reads_head = run_dir + "/" + READS;
  String graph_head = sub_dir + "/" + GRAPH;
  
  String KS = ToString(K);
  
  // Read in read  pairing data.
  cout << Date( ) << ": Reading in files..." << endl;
  
  PairsManager pairs( reads_head + ".pairs");
  
  // Find ploidy.
  const int ploidy = FirstLineOfFile( run_dir + "/ploidy" ).Int( );
  
  vec<int> rlen; // in kmers!
  vec<unipath_id_t> to_rc;
  vec<nkmers_t> ulen;
  int nuni;
  digraph adj_graph;
  
  { // SCOPING TO SAVE MEMORY
    
    // Read in read paths.  We only need the reads' lengths (in kmers).
    vecKmerPath paths(reads_head + ".paths.k" + KS);
    rlen.resize( paths.size() );
    for ( size_t i = 0; i < paths.size(); i++ )
      rlen[i] = paths[i].KmerCount();
    
    // Read in the unipaths.  We only need to know their lengths and their
    // rc mapping.
    vec<tagged_rpint> unipathsdb;
    BinaryRead3( reads_head + ".unipathsdb.k" + KS, unipathsdb );
    vecKmerPath unipaths(reads_head + ".unipaths.k" + KS );
    nuni = unipaths.size( );
    UnipathInvolution( unipaths, unipathsdb, to_rc );
    ulen.resize(nuni);
    for ( int i = 0; i < nuni; i++ )
      ulen[i] = unipaths[i].KmerCount( );
    
    // If we will need the unipath adjacency graph later, we fill it in now,
    // while these data structures are still loaded in memory.
    if ( USE_ADJACENCIES_AS_LINKS ) {
      vecKmerPath pathsrc( reads_head + ".paths_rc.k" + ToString(K) );
      BREAD2( reads_head + ".pathsdb.k" + KS, vec<tagged_rpint>, pathsdb );
      BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths,
				  unipathsdb, adj_graph );
    }
  }
  
  
  
  // Define the normal unipaths.
  // "Normal" here means having a copy number not more than twice the ploidy.
  vec<Bool> normal(nuni, False);
  vec<copy_num_t> predicted_copyno(nuni, -1);
  { // SCOPING TO SAVE MEMORY
    VecPdfEntryVec cp;
    cp.ReadAll( reads_head + ".unipaths.predicted_count.k" + KS );
    ForceAssertEq( (int)cp.size(), nuni );
    
    for ( int i = 0; i < nuni; i++ ) {
      GetMostLikelyValue( predicted_copyno[i], cp[i] );
      
      if ( predicted_copyno[i] > NORMAL_MULT * ploidy ) continue;
      normal[i] = True;
    }
  }

  // Load the read placements on unipaths.
  
  String unilocs_file = reads_head + ".unilocs." + KS;
  BREAD2( unilocs_file, vec<ReadLocationLG>, ulocs );
  VecULongVec ulocs_indexr;
  ulocs_indexr.ReadAll( unilocs_file + ".indexr" );
  
  
  // Sanity checks.  If these asserts fail, the input datasets are mismatched.
  {
    size_t size_sum = 0;
    for ( size_t i = 0; i < ulocs_indexr.size(); i++ )
      size_sum += ulocs_indexr[i].size();
    ForceAssertEq( rlen.size(), pairs.nReads() );
    ForceAssertEq( rlen.size(), ulocs_indexr.size() );
    ForceAssertEq( ulocs.size(), size_sum );
  }
  
  
  /* Build the 2 unipath link graphs
   * Graph 1 ("unipathgraph_seeds"): Used in seed selection.
   * Graph 2 ("unipathgraph_cloud"): Used in cloud creation in local assemblies.
   *
   * The only difference is that, when TRANSITIVE_FILL_IN=True, Graph 1 contains
   * transitive edges (made with FillInTransitiveEdges) while Graph 2 does not.
   *
   */
  
  // Make the graph
  digraphE<fsepdev> LG;
  cout << Date( ) << ": Building unipath link graph" << endl;
  BuildUnipathLinkGraph( K, ploidy, pairs, ulocs, ulocs_indexr, normal,
			 ulen, rlen, to_rc, MIN_EDGE_MULTIPLICITY,
			 predicted_copyno, LG, VERBOSE );
  
  if (PRINT_LINKS)
    for ( int v = 0; v < LG.N( ); v++ ) {
      for ( int i = 0; i < LG.From(v).isize( ); i++ ) {
	int w = LG.From(v)[i];
	fsepdev s = LG.EdgeObjectByIndexFrom( v, i );
	cout << v << " --- " << s.Sep( ) << " +/- " << s.Dev( )
	     << " --> " << w << "\n";
      }
    }
     
  // Add links implied by the unipath graph.
  
  const int max_cell_size = 6;
  if ( USE_ADJACENCIES_AS_LINKS ) {
    cout << Date( ) << ": Adding links that are implied by the unipath "
	 << "adjacency graph" << endl;
    
    vec< vec<int> > cells;
    cout << Date( ) << ": Begin FindCells" << endl;
    FindCells( adj_graph, max_cell_size, cells );
    cout << Date( ) << ": End FindCells" << endl;
    
    for ( int i = 0; i < cells.isize( ); i++ ){    
      vec< vec<int> > allpaths;
      int v = cells[i].front( ), w = cells[i].back( );
      if ( LG.HasEdge( v, w ) ) continue;
      
      // Find all paths from u to v.  This used to be done by
      //      adj_graph.AllPaths( v, w, x );
      // which as it turns out is ridiculously slow by virtue of how
      // AllPaths works, which is to first find all vertices that are
      // upstream of v.  This code is different because we know that
      // every path starting at v must reach w.
      
      set< pair< vec<int>, set<int> > > partials;
      vec<int> just;
      set<int> justs;
      just.push_back(v), justs.insert(v);
      partials.insert( make_pair( just, justs ) );
      set< vec<int> > pathss;
      while( !partials.empty( ) ){    
	pair< vec<int>, set<int> > p = *partials.begin( );
	partials.erase( partials.begin( ) );
	int x = p.first.back( );
	Bool to_process = True;
	if ( x == w ) {    
	  pathss.insert( p.first );
	  to_process = False;    
	}
	if (to_process){    
	  for ( int i = 0; i < adj_graph.From(x).isize( ); i++ ){    
	    int y = adj_graph.From(x)[i];
	    if ( Member( p.second, y ) ) continue;
	    pair< vec<int>, set<int> > q = p;
	    q.first.push_back(y);
	    q.second.insert(y);
	    partials.insert(q);    
	  }    
	}    
      }
      allpaths.reserve( pathss.size( ) );
      for ( set< vec<int> >::iterator i = pathss.begin( ); 
	    i != pathss.end( ); ++i ){    
	allpaths.push_back(*i);    
      }
      
      if ( allpaths.empty( ) ) continue;
      vec<int> dist;
      for ( int i = 0; i < allpaths.isize( ); i++ ){    
	int d = 0;
	for ( int j = 0; j < allpaths[i].isize( ) - 1; j++ )
	  d += rlen[ allpaths[i][j] ];
	dist.push_back(d);    
      }
      int m = Min(dist), M = Max(dist);
      Sort(dist);
      LG.AddEdge( v, w, fsepdev( (m+M)/2, (M-m)/2 ) );    
    }    
  }
  
  int n_links = LG.EdgeObjectCount( );
  // Write the graph object (which does not yet have transitive edges) to unipathgraph_cloud
  // Create the directory for it, if necessary
  
  if (WRITE) {
    Mkdir777( run_dir + "/ASSEMBLIES" );
    Mkdir777( sub_dir );
    BinaryWrite( graph_head + ".cloud.k" + KS, LG );
  }
  
  if (TRANSITIVE_FILL_IN) {
    cout << Date( ) << ": Filling in link graph" << endl;
    
    FillInTransitiveEdges( LG, TRANSITIVE_RADIUS, 1000.0, 10.0,
			   predicted_copyno, ploidy, ulen,
			   MIN_UNIPATH_FOR_TRANSITIVE_JOIN );
  }
  int n_links_transitive = LG.EdgeObjectCount( ) - n_links;
  
  // Write the graph object (now with transitive edges) to unipathgraph_seeds
  if (WRITE)
    BinaryWrite( graph_head + ".seeds.k" + KS, LG );
  
  
  if ( TRANSITIVE_FILL_IN )
    cout << "\tGraph has " << n_links << " links, plus " << n_links_transitive << " transitive edges added." << endl;
  else cout << "\tGraph has " << n_links << " links." << endl;
  
  cout << Date( ) << ": Done!" << endl;
}
