///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/BuildScaffoldGraph.h"
#include "paths/reporting/CLinkBundle.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * LogEvent
 */
void LogEvent( ostream &out,
	       int v1, int v2,
	       int sep, int sd,
	       int weight, float score )
{
  int id1 = ( v1 % 2 == 0 ? v1 / 2 : ( v1 - 1 ) / 2 );
  int id2 = ( v2 % 2 == 0 ? v2 / 2 : ( v2 - 1 ) / 2 );
  bool rc1 = ( v1 % 2 == 1 );
  bool rc2 = ( v2 % 2 == 1 );

  out << "link s_" << id1 << ( rc1 ? "[-]" : "[+]" )
      << "  ->  s_" << id2 << ( rc2 ? "[-]" : "[+]" )
      << "  <" << sep << " +/- " << sd << ">"
      << "  w=" << weight
      << "  score=" << ToString( score, 2 )
      << "\n";
}

/**
 * BuildScaffoldGraph
 */
void BuildScaffoldGraph( const PairsManager &pairs,
			 const vec<superb> &supers,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 digraphE<sepdev> &graph,
			 digraphE<CLinkBundle> *bgraph,
			 vec< pair<int,int> > *black_list,
			 ostream *log,
			 float MAX_SCORE,
			 int MIN_LINKS,
			 float MIN_SCORE,
			 int RATIO_MIN_LINKS,
			 double RATIO_TO_SECOND,
			 int MAX_OVERLAP,
			 int LOW_WEIGHT,
			 double slop )
{
  // Constants.
  const double multiplier = double( RATIO_MIN_LINKS - 1 ) * RATIO_TO_SECOND;
  
  // The log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  out << Date( ) << ": starting BuildScaffoldGraph\n" << endl;

  // Utility maps (only used if black_list is provided).
  vec< vec<int> > nolinks( supers.size( ) );
  if ( black_list ) {
    for (int ii=0; ii<black_list->isize( ); ii++) {
      int s1 = (*black_list)[ii].first;
      int s2 = (*black_list)[ii].second;
      nolinks[s1].push_back( s2 );
      nolinks[s2].push_back( s1 );
    }
  }
  
  // The linker.
  CSuperLinks linker( &pairs, &supers, &aligns, &index );
  int n_supers = supers.size( );
  
  // Core data to generate a digraphE.
  vec< vec<int> > from( 2 * n_supers );
  vec< vec<int> > to( 2 * n_supers );
  vec< vec<int> > from_edge_obj( 2 * n_supers );
  vec< vec<int> > to_edge_obj( 2 * n_supers );
  vec<sepdev> edges;
  vec<CLinkBundle> bedges;
  
  // Loop over all supers.
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    vec<COffset> all_offsets;
    linker.AllLinks( super_id, all_offsets, slop );

    // Verbose log (print all links).
    if ( log ) 
      for (int ii=0; ii<all_offsets.isize( ); ii++)
      	all_offsets[ii].Print( out, &pairs );
    
    // Avoid duplication. Pool together super2s with different orientations.
    vec< vec<COffset> > offsets;
    for (size_t ii=0; ii<all_offsets.size( ); ii++) {
      int s2 = all_offsets[ii].Super2( );
      ForceAssertEq( all_offsets[ii].Super1( ), super_id );
      ForceAssertNe( super_id, s2 );
      if ( super_id > s2 ) continue;
      if ( black_list && ! ( all_offsets[ii].Rc2( ) ) ) {
	bool skip = false;
	for (int jj=0; jj<nolinks[super_id].isize( ); jj++)
	  if ( nolinks[super_id][jj] == s2 ) {
	    skip = true;
	    break;
	  }
	if ( skip ) continue;
      }
      if ( offsets.size( ) < 1 || offsets.back( ).back( ).Super2( ) != s2 )
	offsets.resize( offsets.size( ) + 1 );
      offsets[ offsets.size( ) - 1 ].push_back( all_offsets[ii] );
    }

    // For every pair of linked supers, decided if there is a winner cluster.
    for (size_t mate_id=0; mate_id<offsets.size( ); mate_id++) {
      ForceAssertLe( offsets[mate_id].isize( ), 2 );

      // Container for combined_scores, sizes and signed cluster ids.
      vec< triple<double,int,int> > organizer;   // (score, weight, signed_id)
      for (int orient_id=0; orient_id<offsets[mate_id].isize( ); orient_id++) {
	const COffset &offset = offsets[mate_id][orient_id];
	for (size_t ii=0; ii<offset.NClusters( ); ii++) {
	  int nlinks = offset.NLinks( ii );
	  float score = offset.Score( ii );
	  if ( nlinks < 2 ) continue;
	  if ( nlinks < MIN_LINKS && score < MIN_SCORE ) continue;
	  if ( score > MAX_SCORE ) continue;
	  int signed_id = orient_id == 0 ? (int)ii : - 1 - (int)ii;
	  double comb_score = offset.CombinedScore( ii );
	  triple<double,int,int> newtriple( comb_score, nlinks, signed_id );
	  organizer.push_back( newtriple );
	}    
      }
      
      // Sort organizers (heavier clusters go first).
      sort( organizer.rbegin( ), organizer.rend( ) );
      if ( organizer.size( ) < 1 ) continue;
      
      // There is clear winner.
      bool is_good = false;
      if ( organizer.size( ) == 1 ) is_good = true;
      else if ( organizer.size( ) > 1 ) {
	double size0 = organizer[0].second;
	double size1 = organizer[1].second;
	if ( size0 > multiplier * size1 ) is_good = true;
      }
      if ( ! is_good ) continue;

      // Use winning cluster to link supers.
      int wcluster_id = organizer[0].third;
      const COffset &winner = offsets[mate_id][ wcluster_id < 0 ? 1 : 0 ];
      if ( wcluster_id < 0 ) wcluster_id = - 1 - wcluster_id;
      
      const int weight = winner.NLinks( wcluster_id );
      const float score = winner.Score( wcluster_id );
      const normal_distribution ndoff = winner.Offset( wcluster_id );
      const pair<double,double> win1 = winner.SpreadWin1( wcluster_id );
      const pair<double,double> win2 = winner.SpreadWin2( wcluster_id );
      if ( ndoff.mu_ >= 0 ) {
	int sep = (int)ndoff.mu_ - winner.Slen1( );
	int sd = (int)ndoff.sigma_;
	
	int v1 = 2 * winner.Super1( );
	int v2 = 2 * winner.Super2( ) + ( winner.Rc2( ) ? 1 : 0 );
	from[v1].push_back( v2 );
	to[v2].push_back( v1 );
	from_edge_obj[v1].push_back( edges.isize( ) );
	to_edge_obj[v2].push_back( edges.isize( ) );
	edges.push_back( sepdev( sep, sd ) );
	bedges.push_back( CLinkBundle( sep, sd, weight, score, win1, win2 ) );
	ForceAssertNe( v1, v2 );
	// LogEvent( out, v1, v2, sep, sd, weight, score );
	
	v1 = 2 * winner.Super2( ) + ( winner.Rc2( ) ? 0 : 1 );
	v2 = 2 * winner.Super1( ) + 1;
	from[v1].push_back( v2 );
	to[v2].push_back( v1 );
	from_edge_obj[v1].push_back( edges.isize( ) );
	to_edge_obj[v2].push_back( edges.isize( ) );
	edges.push_back( sepdev( sep, sd ) );
	bedges.push_back( CLinkBundle( sep, sd, weight, score, win2, win1 ) );
	ForceAssertNe( v1, v2 );
	// LogEvent( out, v1, v2, sep, sd, weight, score );
      }
      else {
	int sep = - (int)ndoff.mu_ - winner.Slen2( );
	int sd = (int)ndoff.sigma_;
	
	int v1 = 2 * winner.Super2( ) + ( winner.Rc2( ) ? 1 : 0 );
	int v2 = 2 * winner.Super1( );
	from[v1].push_back( v2 );
	to[v2].push_back( v1 );
	from_edge_obj[v1].push_back( edges.isize( ) );
	to_edge_obj[v2].push_back( edges.isize( ) );
	edges.push_back( sepdev( sep, sd ) );
	bedges.push_back( CLinkBundle( sep, sd, weight, score, win2, win1 ) );
	ForceAssertNe( v1, v2 );
	// LogEvent( out, v1, v2, sep, sd, weight, score );

	v1 = 2 * winner.Super1( ) + 1;
	v2 = 2 * winner.Super2( ) + ( winner.Rc2( ) ? 0 : 1 );
	from[v1].push_back( v2 );
	to[v2].push_back( v1 );
	from_edge_obj[v1].push_back( edges.isize( ) );
	to_edge_obj[v2].push_back( edges.isize( ) );
	edges.push_back( sepdev( sep, sd ) );
	bedges.push_back( CLinkBundle( sep, sd, weight, score, win1, win2 ) );
	ForceAssertNe( v1, v2 );
	// LogEvent( out, v1, v2, sep, sd, weight, score );
      }
    }
  }

  // Vectors from and to need to be sorted.
  for (int ii=0; ii<2*n_supers; ii++) {
    sort( from[ii].begin( ), from[ii].end( ) );
    sort( to[ii].begin( ), to[ii].end( ) );
  }
  
  // Generate digraph.
  out << Date( ) << ": generating graph" << endl;
  graph.Initialize( from, to, edges, to_edge_obj, from_edge_obj );
  if ( bgraph )
    bgraph->Initialize( from, to, bedges, to_edge_obj, from_edge_obj );

  // Most of the edges with large overlap and low weight are innies.
  vec<int> to_left;
  bgraph->ToLeft( to_left );
  int n_edges = bgraph->EdgeObjectCount( );
  
  out << Date( ) << ": " << n_edges << " edges found" << endl;
  vec<int> to_delete;
  for (int ii=0; ii<n_edges; ii++) {
    int sep = bgraph->EdgeObject( ii ).Sep( );
    int weight = bgraph->EdgeObject( ii ).Weight( );
    if ( sep < - MAX_OVERLAP && weight < LOW_WEIGHT )
      to_delete.push_back( ii );
  }
  
  out << Date( ) << ": cutting "
      << to_delete.size( ) << " edges (possibly innies)"
      << endl; 
  bgraph->DeleteEdges( to_delete, to_left );
  bgraph->RemoveDeadEdgeObjects( );
  graph.DeleteEdges( to_delete, to_left );
  graph.RemoveDeadEdgeObjects( );
  
  // Log final graph.
  if ( log ) {
    out << "\nOUTPUT GRAPH:\n" << endl;

    vec<int> to_left;
    vec<int> to_right;
    bgraph->ToLeft( to_left );
    bgraph->ToRight( to_right );
    
    for (int edge_id=0; edge_id<bgraph->EdgeObjectCount( ); edge_id++) {
      int v1 = to_left[edge_id];
      int s1 = v1 / 2;
      int len1 = supers[s1].TrueLength( );;
      bool rc1 = ( 0 != v1 % 2 );
      
      int v2 = to_right[edge_id];
      int s2 = v2 / 2;
      int len2 = supers[s2].TrueLength( );;
      bool rc2 = ( 0 != v2 % 2 );
      
      out << "s" << s1 << ( rc1 ? "[-]" : "[+]" ) << " (" << len1 << " bp)"
	  << "  to  "
	  << "s" << s2 << ( rc2 ? "[-]" : "[+]" ) << " (" << len2 << " bp)"
	  << "  :  " << bgraph->EdgeObject( edge_id ).AsString( true )
	  << "\n";
    }

    out << endl;
  }

  // Done.
  out << Date( ) << ": BuildScaffoldGraph done" << endl;
  
}
