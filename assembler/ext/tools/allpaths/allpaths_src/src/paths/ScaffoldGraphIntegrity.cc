///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Superb.h"
#include "Vec.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/ScaffoldGraphIntegrity.h"

/**
 * LogSupers
 */
String LogSupers( int id1, int id2 )
{
  int s1 = id1 / 2;
  int s2 = id2 / 2;
  bool rc1 = ( 0 != id1 % 2 );
  bool rc2 = ( 0 != id2 % 2 );

  return
    "s" + ToString( s1 ) + ( rc1 ? "[-]" : "[+]" )
    + " to s" + ToString( s2 ) + ( rc2 ? "[-]" : "[+]" );
}

/**
 * ScaffoldGraphIntegrity
 */
bool ScaffoldGraphIntegrity( const digraphE<CLinkBundle> &graph, ostream *log )
{
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;

  // Total counter for problems found.
  int n_errors = 0;

  // Maps and such.
  int n_edges = graph.EdgeObjectCount( );
  int n_vertices = graph.N( );
  vec<int> to_left;
  vec<int> to_right;
  graph.ToLeft( to_left );
  graph.ToRight( to_right );

  // Identify twin edges, and check they match.
  vec<int> twin_of( n_edges, -1 );

  for (int idv=0; idv<n_vertices; idv++) {
    vec<int> all_idws = graph.From( idv );
    for (int jj=0; jj<all_idws.isize( ); jj++) {
      int idw = all_idws[jj];
      vec<int> edges = graph.EdgesBetween( idv, idw );
      if ( edges.size( ) != 1 ) {
	n_errors++;
	out << edges.size( ) << " edge(s) found from "
	    << idv << " and "
	    << idw << " ((it should be 1)\n";
	continue;
      }
      int rc_idv = ( 0 == idv % 2 ) ? idv + 1 : idv - 1;
      int rc_idw = ( 0 == idw % 2 ) ? idw + 1 : idw - 1;
      vec<int> rc_edges = graph.EdgesBetween( rc_idw, rc_idv );
      if ( rc_edges.size( ) != 1 ) {
	n_errors++;
	out << rc_edges.size( ) << " edge(s) found from "
	    << rc_idw << " and "
	    << rc_idv << " ((it should be 1)\n";
	continue;
      }
      twin_of[ edges[0] ] = rc_edges[0];
      twin_of[ rc_edges[0] ] = edges[0];
      const CLinkBundle &bundle1 = graph.EdgeObject( edges[0] );
      const CLinkBundle &bundle2 = graph.EdgeObject( rc_edges[0] );
      if ( bundle1.Sep( ) != bundle2.Sep( ) ||
	   bundle1.Dev( ) != bundle2.Dev( ) ||
	   bundle1.Weight( ) != bundle2.Weight( ) ||
	   bundle1.Score( ) != bundle2.Score( ) ) {
	n_errors++;
	out << "Broken symmetry: "
	    << LogSupers( idv, idw )
	    << " [" << bundle1.AsString( true ) << "], but "
	    << LogSupers( rc_idw, rc_idv )
	    << " [" << bundle2.AsString( true ) << "]\n";
      }
    }
  }
  
  // Check all edges are accounted for.
  for (int ii=0; ii<twin_of.isize( ); ii++) {
    if ( twin_of[ii] < 0 ) {
      n_errors++;
      out << "Edge " << ii << " does not have valid twin\n";
    }
  }


  // Log total tally and leave.
  if ( n_errors > 0 ) {
    out << n_errors << " ERRORS FOUND! ScaffoldGraphIntegrity failed" << endl;
    return false;
  }
  out << "ScaffoldGraphIntegrity run: no errors found" << endl;
  return true;
  
}

