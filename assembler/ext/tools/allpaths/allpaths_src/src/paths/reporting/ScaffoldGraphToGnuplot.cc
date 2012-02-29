/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Superb.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "paths/reporting/CLinkBundle.h"

/**
 * ScaffoldGraphToGnuplot
 *
 * Load a given scaffold graph, and print all its edges in a gnuplot
 * friend format. It needs in input HEAD.super and HEAD.graph.
 */
int main( int argc, char *argv[] )
{
  BeginCommandArguments;
  CommandArgument_String( HEAD );
  EndCommandArguments;

  // Dir and file names.
  String supers_file = HEAD + ".superb";
  String graph_file = HEAD + ".graph";
  String out_file = HEAD + ".graph.txt";

  ofstream out( out_file.c_str( ) );
  PrintCommandPretty( out );
  
  cout << "Sending output to " << out_file << "\n" << endl;

  // Load.
  out << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  out << Date( ) << ": loading graph\n" << endl;
  digraphE<CLinkBundle> graph;
  int fd = OpenForRead( graph_file );
  BinaryRead( fd, graph );
  
  // Build indices.
  vec<int> to_left;
  vec<int> to_right;
  graph.ToLeft( to_left );
  graph.ToRight( to_right );
  
  // Legend.
  out << "# "
      << "id1 "
      << "size1 "
      << "id2 "
      << "size2 "
      << "sep "
      << "dev "
      << "weight "
      << "score "
      << "begin1 "
      << "end1 "
      << "begin2 "
      << "end2\n";
  
  // Loop over all edges.
  int n_edges = graph.EdgeObjectCount( );
  for (int edge_id=0; edge_id<n_edges; edge_id++) {
    const CLinkBundle &edge = graph.EdgeObject( edge_id );
    ForceAssertNe( to_left[edge_id], to_right[edge_id] );
    int id1 = to_left[edge_id] / 2;
    int id2 = to_right[edge_id] / 2;
    bool rc1 = 1 == to_left[edge_id] % 2;
    bool rc2 = 1 == to_right[edge_id] % 2;
    int len1 = ( rc1 ? -1 : +1 ) * supers[id1].TrueLength( );
    int len2 = ( rc2 ? -1 : +1 ) * supers[id2].TrueLength( );
    pair<double,double> win1 = edge.Win1( );
    pair<double,double> win2 = edge.Win2( );
    
    out << id1 << " "
	<< len1 << " "
	<< id2 << " "
	<< len2 << " "
	<< edge.Sep( ) << " "
	<< edge.Dev( ) << " "
	<< edge.Weight( ) << " "
	<< edge.Score( ) << " "
	<< win1.first << " "
	<< win1.second << " "
	<< win2.first << " "
	<< win2.second << "\n";
  }
  out << endl;
  
  // Done.
  cout << Date( ) << ": done" << endl;
  out << Date( ) << ": done" << endl;
  out.close( );

}

