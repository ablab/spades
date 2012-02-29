/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/Sepdev.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/MapAligns.h"
#include "paths/MakeScaffoldsScoredBest.h"
#include "paths/ScaffoldGraphIntegrity.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/SuckScaffolds.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * PrintInfo
 *
 * Helper function to print core info for an edge between two supers.
 */
void PrintInfo( ostream &out,
		const int v1,
		const int v2,
		const CLinkBundle &edge,
		const vec<superb> &supers )
{
  int s1 = v1 / 2;
  int len1 = supers[s1].TrueLength( );;
  bool rc1 = ( 0 != v1 % 2 );

  int s2 = v2 / 2;
  int len2 = supers[s2].TrueLength( );;
  bool rc2 = ( 0 != v2 % 2 );
  
  out << "s" << s1 << ( rc1 ? "[-]" : "[+]" ) << " (" << len1 << " bp)"
      << "  to  s" << s2 << ( rc2 ? "[-]" : "[+]" ) << " (" << len2 << " bp)"
      << "  :  " << edge.AsString( true )
      << "\n";
  
}

/**
 * AddToCandidates
 *
 * Append edge ids to the lists of candidate keepers and candidate
 * deleters.
 */
void AddToCandidates( const digraphE<CLinkBundle> &graph,
		      const vec<int> &select,
		      vec<int> &candidate_keepers,
		      vec<int> &candidate_deleters )
{
  for (size_t ii=0; ii<select.size( ); ii++) {
    int edge_to = -1;
    if ( ii > 0 ) {
      int v = select[ii-1];
      int w = select[ii];
      vec<int> edges = graph.EdgesBetween( v, w );
      edge_to = edges[0];
    }
    
    int edge_from = -1;
    if ( ii < select.size( ) - 1 ) {
      int v = select[ii];
      int w = select[ii+1];
      vec<int> edges = graph.EdgesBetween( v, w );
      edge_from = edges[0];
      candidate_keepers.push_back( edge_from );
    }
    
    if ( ii > 0 ) {
      int v = select[ii];
      vec<int> all_tos = graph.ToEdgeObj( v );
      for (int ii=0; ii<all_tos.isize( ); ii++)
	if ( all_tos[ii] != edge_to )
	  candidate_deleters.push_back( all_tos[ii] );
    }
    
    if ( ii < select.size( ) - 1 ) {
      int v = select[ii];
      vec<int> all_froms = graph.FromEdgeObj( v );
      for (int ii=0; ii<all_froms.isize( ); ii++)
	if ( all_froms[ii] != edge_from )
	  candidate_deleters.push_back( all_froms[ii] );
    }
  }
  
}

/**
 * FindPaths
 *
 * Find paths between the two vertices, with some constraints (to
 * avoid run-time explosions on tangled parts of the scaffold
 * graph). Just a wrapper around digraphE<E>::AllPathsLengthRange.
 *
 * It assumes that v and w are connected by a direct edge, and the
 * output will always contain at least the direct path v1 to v2.
 */
void FindPaths( const int v,
		const int w,
		const digraphE<int> &igraph,
		const vec<int> &to_right,
		vec< vec<int> > &paths )
{
  // Clear.
  paths.clear( );

  // Range for lengths.
  const int L1 = 0;
  const int L2 = 8;

  // Call AllPathsLengthRange.
  vec< vec<int> > edge_paths;
  bool ok = igraph.AllPathsLengthRange( v, w, L1, L2, to_right, edge_paths );

  if ( ok ) {   // Turn paths of edges into paths of vertices
    for (int ii=0; ii<edge_paths.isize( ); ii++) {
      vec<int> vpath( 1, v );
      for (int jj=0; jj<edge_paths[ii].isize( ); jj++)
	vpath.push_back( to_right[ edge_paths[ii][jj] ] );
      paths.push_back( vpath );
    }
  }
  else {   // No paths found, default to direct link
    if ( paths.size( ) < 1 ) {
      vec<int> direct;
      direct.push_back( v );
      direct.push_back( w );
      paths.push_back( direct );
    }
  }

}

/**
 * MakeScaffoldsScoredBest
 */
void MakeScaffoldsScoredBest( const PairsManager &pairs,
			      vec<fastavector> &contigs,
			      vec<superb> &supers,
			      vec<alignlet> &aligns,
			      vec<int> &index,
			      ostream &out,
			      int *MIN_LINKS,
			      int *MIN_SEP,
			      bool SUCK_SCAFFOLDS,
			      bool VERBOSE )
{
  // HEURISTICS - used in the triangularization process.
  const double MAX_STRETCH = 10.0;

  // Build graph.
  if ( VERBOSE ) {
    out << Date( ) << ": starting MakeScaffoldsScoredBest\n"
	<< Date( ) << ": building graph" << endl;
  }
  digraphE<sepdev> unused;
  digraphE<CLinkBundle> graph;
  BuildScaffoldGraph( pairs, supers, aligns, index, unused, &graph );
  ForceAssert( ScaffoldGraphIntegrity( graph ) );
  int n_edges_orig = graph.EdgeObjectCount( );

  // Delete edges from graph.
  if ( n_edges_orig > 0 ) {
    vec<int> to_delete;
    to_delete.reserve( n_edges_orig );
    
    // Delete test: remove edges with low weight.
    if ( MIN_LINKS ) {
      for (int ii=0; ii<n_edges_orig; ii++)
	if ( graph.EdgeObject( ii ).Weight( ) < *MIN_LINKS ) 
	  to_delete.push_back( ii );
    }
    
    // Delete test: remove edges with large overlap.
    if ( MIN_SEP ) {
      for (int ii=0; ii<n_edges_orig; ii++)
	if ( graph.EdgeObject( ii ).Sep( ) < *MIN_SEP )
	  to_delete.push_back( ii );
    }
    
    // Actually remove edges from graph (and regenerate to_left!)
    vec<int> to_left;
    graph.ToLeft( to_left );
    graph.DeleteEdges( to_delete, to_left );
    graph.RemoveDeadEdgeObjects( );
    graph.ToLeft( to_left );
    
    int n_edges_after = graph.EdgeObjectCount( );
    int n_removed = n_edges_orig - n_edges_after;
    double ratio = 100.0 * SafeQuotient( n_removed, n_edges_orig );
    if ( VERBOSE )
      out << "removing edges:\n"
	  << n_edges_orig << " edges in input, "
	  << n_edges_after << " in output ("
	  << n_removed << " edges removed, correspondig to "
	  << ToString( ratio, 2 ) << "%)\n";
  }
  
  // No edges left in graph. Run the original MakeScaffolds and return.
  int n_edges = graph.EdgeObjectCount( );
  if ( n_edges < 1 ) return;
  
  if ( VERBOSE ) 
    out << Date( ) << ": " << n_edges << " edges left in graph" << endl;
  
  // Maps to nagavigate the initial graph.
  int n_vertices = graph.N( );
  vec<int> to_left;
  vec<int> to_right;
  graph.ToLeft( to_left );
  graph.ToRight( to_right );

  // Sort bundles.
  vec<const CLinkBundle*> bundles;
  vec<int> bids;
  bundles.reserve( n_edges );
  bids.reserve( n_edges );
  for (int ii=0; ii<n_edges; ii++) {
    bundles.push_back( &( graph.EdgeObject( ii ) ) );
    bids.push_back( ii );
  }
  pCLinkBundle_order_best sorter;
  SortSync( bundles, bids, sorter );
  
  if ( VERBOSE ) {
    out << "best bundles found:\n";
    for (int ii=0; ii<bundles.isize( ); ii++) {
      int idv = to_left[ bids[ii] ];
      int idw = to_right[ bids[ii] ];
      const CLinkBundle *pbundle = bundles[ii];
      PrintInfo( out, idv, idw, *pbundle, supers );
    }
    if ( bundles.isize( ) < 1 ) out << "no bundles found\n";
  }
  
  
#ifdef SKIP_THIS // SANTEMP
  String out_base = "/wga/scr2/sante/output/temp";

  String graph_file = out_base + ".graph";
  int fd = OpenForWrite( graph_file );
  BinaryWrite( fd, graph );
  
  String supers_file = out_base + ".superb";
  WriteSuperbs( supers_file, supers );
  
  String supers_txt_file = out_base + ".superb.txt";
  String theCommand = "ShowSupersBrief SUPERBS=" + supers_file;
  String theLog = " >& " + supers_txt_file;
  SystemSucceed( theCommand + theLog );

  String links_file = out_base + ".links";
  ofstream lout( links_file.c_str( ) );
  for (int ii=0; ii<bundles.isize( ); ii++) {
    int idv = to_left[ bids[ii] ];
    int idw = to_right[ bids[ii] ];
    const CLinkBundle *pbundle = bundles[ii];
    lout << ii << "\t";
    PrintInfo( lout, idv, idw, *pbundle, supers );
  }
  lout.close( );

  ForceAssert( 1 == 0 );
#endif // SANTEMP
  
  
  // Build a parallel graph with E=int (needed for AllPathsLengthRange).
  vec<int> iedges( graph.Edges( ).size( ), 1 );
  digraphE<int> igraph( graph.From( ),  graph.To( ), iedges,
			graph.ToEdgeObj( ), graph.FromEdgeObj( ) );
  
  // Loop over all bundles (starting from the best ones).
  vec<bool> to_keep( n_edges, false );
  vec<bool> to_delete( n_edges, false );
  for (int bpos=0; bpos<bundles.isize( ); bpos++) {
    const int edge_id = bids[bpos];
    const CLinkBundle *p_edge = bundles[bpos];
    const int vtx_v = to_left[edge_id];
    const int vtx_w = to_right[edge_id];
    
    // Paths between v and w.
    vec< vec<int> > paths;
    FindPaths( vtx_v, vtx_w, igraph, to_right, paths );
    
    
#ifdef SKIP_THIS // SANTEMP
    if ( vtx_v == 1318 && vtx_w == 373 ) {
      cout << "BEFORE:\n";
      for (int ii=0; ii<paths.isize( ); ii++) {
	for (int jj=0; jj<paths[ii].isize( ); jj++)
	  cout << paths[ii][jj] << " ";
	cout << "\n";
      }
      cout << endl;
    }
#endif // SANTEMP
    
    
    // Try to merge paths, whenever possible.
    {
      vec<int> plens( paths.size( ) );
      for (int ii=0; ii<paths.isize( ); ii++)
	plens[ii] = paths[ii].isize( );
      ReverseSortSync( plens, paths );
      
      vec<int>::iterator it;
      for (int ii=1; ii<paths.isize( )-1; ii++) {   // -1: keep direct path.
	bool is_embedded = false;
	for (int jj=0; jj<ii; jj++) {
	  bool is_missing = false;
	  for (int kk=1; kk<paths[ii].isize( )-1; kk++) {
	    it = find( paths[jj].begin( ), paths[jj].end( ), paths[ii][kk] );
	    if ( it == paths[jj].end( ) ) {
	      is_missing = true;
	      break;
	    }
	  }
	  if ( ! is_missing ) {
	    is_embedded = true;
	    break;
	  }
	}
	if ( is_embedded ) {
	  paths.erase( paths.begin( ) + ii );
	  ii--;
	}
      }
    }
    
    
#ifdef SKIP_THIS // SANTEMP
    if ( vtx_v == 1318 && vtx_w == 373 ) {
      cout << "AFTER:\n";
      for (int ii=0; ii<paths.isize( ); ii++) {
	for (int jj=0; jj<paths[ii].isize( ); jj++) 
	  cout << paths[ii][jj] << " ";
	cout << "\n";
      }
      cout << endl;
      ForceAssert( 1 == 0 );
    }
#endif // SANTEMP
    
    
    // The final path (for now skip complicate cases).
    const vec<int> *select = 0;
    if ( paths.size( ) != 2 ) {
      for (size_t ii=0; ii<paths.size( ); ii++) {
	if ( paths[ii].size( ) == 2 ) {
	  select = &( paths[ii] );
	  break;
	}
      }
    }
    else {
      const vec<int> &direct = paths[0].size( ) == 2 ? paths[0] : paths[1];
      const vec<int> &full = paths[0].size( ) == 2 ? paths[1] : paths[0];
      
      // Is the longest path consistent with the direct (single edge) one?
      int direct_sep = p_edge->Sep( );
      int direct_dev = p_edge->Dev( );
      int direct_rad = ( MAX_STRETCH * (double)direct_dev );
      ho_interval directw( direct_sep - direct_rad, direct_sep + direct_rad );

      int full_sep = 0;
      int full_dev_square = 0;
      vec<CLinkBundle> bundles;
      for (int ii=0; ii<full.isize( )-1; ii++) {
	bundles = graph.EdgeObjectsBetween( full[ii], full[ii+1] );
	int n_bundles = bundles.size( );
	ForceAssertEq( n_bundles, 1 );
	full_sep += bundles[0].Sep( );
	full_dev_square += bundles[0].Dev( ) * bundles[0].Dev( );
	if ( ii > 0 ) full_sep += supers[ full[ii] / 2 ].TrueLength( );
      }
      int full_dev = sqrt( (double)full_dev_square );
      int full_rad = ( MAX_STRETCH * (double)full_dev );
      ho_interval fullw( full_sep - full_rad, full_sep + full_rad );
      

#ifdef SKIP_THIS // SANTEMP
      if ( vtx_v == 1318 && vtx_w == 373 ) {
	cout << "IN ++++++++++++++++++++++++++++\n";
	cout << "direct: <" << direct_sep << " +/- " << direct_dev << ">\n";
	cout << "full: <" << full_sep << " +/- " << full_dev << ">\n";
	cout << ".....................\n";
	cout << "direct ho: " << directw << "\n";
	cout << "full ho:   " << fullw << "\n";
	cout << "OUT +++++++++++++++++++++++++++\n";
	ForceAssert( 1 == 0 );
      }
#endif // SANTEMP
      
      
      bool is_consistent = ( Overlap( directw, fullw ) > 0 );
      select = is_consistent ? &full : &direct;
    }
    
    if ( ! select ) continue;

    // Candidate keepers sand deleters.
    vec<int> candidate_keepers;
    vec<int> candidate_deleters;
    AddToCandidates( graph, *select, candidate_keepers, candidate_deleters );
    
    vec<int> select_rc;
    select_rc.reserve( select->size( ) );
    for (int ii=select->isize( )-1; ii>=0; ii--) {
      int id = (*select)[ii];
      if ( id % 2 == 0 ) id += 1;
      else id += -1;
      select_rc.push_back( id );
    }
    AddToCandidates( graph, select_rc, candidate_keepers, candidate_deleters );
    
    // Update to_keep and to_delete (if consistent).
    bool is_consistent = true;
    for (int ii=0; ii<candidate_keepers.isize( ); ii++) {
      if ( to_delete[ candidate_keepers[ii] ] ) {
	is_consistent = false;
	break;
      }
    }
    if ( is_consistent ) {
      for (int ii=0; ii<candidate_deleters.isize( ); ii++) {
	if ( to_keep[ candidate_deleters[ii] ] ) {
	  is_consistent = false;
	  break;
	}
      }
    }
    if ( ! is_consistent ) continue;

    for (int ii=0; ii<candidate_keepers.isize( ); ii++)
      to_keep[ candidate_keepers[ii] ] = true;
    for (int ii=0; ii<candidate_deleters.isize( ); ii++)
      to_delete[ candidate_deleters[ii] ] = true;
    
  }
  
  // Remove edges.
  vec<int> del_ids;
  del_ids.reserve( to_delete.size( ) );
  for (size_t ii=0; ii<to_delete.size( ); ii++)
    if ( to_delete[ii] ) del_ids.push_back( (int)ii );
  graph.DeleteEdges( del_ids, to_left );
  graph.RemoveDeadEdgeObjects( );
  graph.ToLeft( to_left );
  
  // Clean the graph by making sure each vertex has as most one to/from.
  vec<int> xdel;
  for (int ii=0; ii<graph.N( ); ii++) {
    vec<int> all_tos = graph.ToEdgeObj( ii );
    if ( all_tos.size( ) > 1 )
      copy( all_tos.begin( ), all_tos.end( ), back_inserter( xdel ) );
    
    vec<int> all_froms = graph.FromEdgeObj( ii );
    if ( all_froms.size( ) > 1 )
      copy( all_froms.begin( ), all_froms.end( ), back_inserter( xdel ) );
  }
  
  if ( xdel.size( ) > 0 ) {
    sort( xdel.begin( ), xdel.end( ) );
    xdel.erase( unique( xdel.begin( ), xdel.end( ) ), xdel.end( ) );
    if ( VERBOSE )
      out << Date( ) << ": removing " << xdel.size( ) << " extra edges" << endl;
    graph.DeleteEdges( xdel, to_left );
    graph.RemoveDeadEdgeObjects( );
    graph.ToLeft( to_left );
  }
  
  // Find sources and build supers.
  vec<int> sources;
  graph.Sources( sources );
  
  vec<superb> old_supers;   // will be swapped and turn into the old supers
  vec<bool> grabbed( supers.size( ), false );
  vec<bool> reversed( contigs.size( ), false );
  old_supers.reserve( sources.size( ) );
  for (int source_id=0; source_id<sources.isize( ); source_id++) {
    int vid = sources[source_id];

    // Avoid duplication.
    if ( grabbed[ vid / 2 ] ) continue;
    vec<int> vpath( 1, vid );

    // Build new super (tag contigs to rc).
    superb sup;
    while ( 1 ) {
      int original_super_id = vid / 2;
      superb original_super = supers[original_super_id];
      if ( 1 == vid % 2 ) {
	original_super.Reverse( );
	for (int jj=0; jj<original_super.Ntigs( ); jj++) {
	  int cg_id = original_super.Tig( jj );
	  reversed[cg_id] = ! reversed[cg_id];
	}
      }
      if ( sup.Ntigs( ) < 1 ) {
	grabbed[original_super_id] = true;
	sup = original_super;
      }
      else {
	vec<int> tos = graph.To( vid );
	ForceAssertEq( tos.isize( ), 1 );
	vec<int> edge_ids = graph.EdgesBetween( tos[0], vid );
	ForceAssertEq( edge_ids.isize( ), 1 );
	const CLinkBundle &edge = graph.EdgeObject( edge_ids[0] );
	
	sup.AppendSuper( original_super, edge.Sep( ), edge.Dev( ) );
	grabbed[original_super_id] = true;
	vpath.push_back( vid );
      }
      vec<int> froms = graph.From( vid );

      // Iterate.
      if ( froms.isize( ) == 1 ) {
	vid = froms[0];
	continue;
      }

      // This should never happen (still generate super, but log event).
      if ( froms.isize( ) > 1 ) out << "WARNING ERROR: multiple next supers. ";
 
      // Break.
      break;
    }
    
    // Add super.
    old_supers.push_back( sup );
  }
  
  // don't forget about the unconnected supers
  for (int s=0; s < supers.isize( ); s++) {
    if (!grabbed[s])
      old_supers.push_back(supers[s]);
  }


  swap( supers, old_supers );
  
  // Flip contigs and aligns.
  if ( VERBOSE ) out << "merging done, flipping rc-ed contigs" << endl;
  for (int contig_id=0; contig_id<reversed.isize( ); contig_id++)
    if ( reversed[contig_id] )
      contigs[contig_id].ReverseComplement( );

  for (size_t ii=0; ii<aligns.size( ); ii++)
    if ( reversed[ aligns[ii].TargetId( ) ] )
      aligns[ii].Reverse( );
  
  // Suck scaffolds. Warning! aligns are flipped in SuckScaffolds!
  if ( SUCK_SCAFFOLDS ) {
    vec<Bool> rctig(contigs.size(), false);
    int n_sucked = SuckScaffolds( pairs, index, aligns,
				  supers, rctig, VERBOSE ? &out : NULL );
    for (u_int c = 0; c < contigs.size(); ++c)
      if ( rctig[c] ) contigs[c].ReverseComplement();
    if ( VERBOSE ) 
      out << Date( ) << ": " << n_sucked << " scaffolds embedded" << endl;
  }

  // Log summary and leave.
  if ( VERBOSE ) {
    int n_supers = supers.isize( );
    vec<int> slens;
    slens.reserve( supers.size( ) );
    for (int ii=0; ii<supers.isize( ); ii++)
      slens.push_back( supers[ii].TrueLength( ) );
    sort( slens.begin( ), slens.end( ) );
    int n50 = ( slens.size( ) > 0 ) ? N50( slens ) : 0;
    
    out << Date( ) << ": "
	<< old_supers.size( ) << " supers in, "
	<< supers.size( ) << " supers out (N50="
	<< n50 << ")" << endl;
    
    out << "\nStats after ScoredBest:\n" << endl;
    ReportScaffoldsN50( supers, out );
    
    out << Date( ) << ": MakeScaffoldsScoredBest done\n" << endl;
  }
  
}

