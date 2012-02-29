///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "graph/Digraph.h"
#include "feudal/IncrementalWriter.h"
#include "paths/GetNexts.h"
#include "paths/KmerPath.h"
#include "PairsManager.h"

/**
 * SplitUnibases
 *
 * Split the given fastb into overlapping chunks, with some rules:
 *   1. bvec's are assumed to have length >= K
 *   2. consecutive chunks overlap by exactly K bases
 *   3. bvec's are broken into chunks CHUNK_LEN long (or less)
 *
 * Input files:
 *   <IN_HEAD>.unibases.k<K>
 *   <IN_HEAD>.unipaths.k<K>
 *   <IN_HEAD>.unipath_adjgraph.k<K>
 *
 * Output files:
 *   <OUT_HEAD>.SplitUnibases.log
 *   <OUT_HEAD>.fastb
 *   <OUT_HEAD>.paths.k<K>
 *   <OUT_HEAD>.pairs* (empty pairing info)
 *
 * K: consecutive chunks overlap by bases
 * IN_HEAD: full path name for input files
 * OUT_HEAD: full path base name for output files
 * CHUNK_LEN: length of chunks (maximum)
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_Int( K );
  CommandArgument_String( IN_HEAD );
  CommandArgument_String( OUT_HEAD );
  CommandArgument_Bool_OrDefault_Doc( USE_ADJ_GRAPH, True,
    "Use pre-computed adj graph or else try K-1 base overlaps." );
  CommandArgument_Int_OrDefault( CHUNK_LEN, 180 );
  EndCommandArguments;

  // File names.
  String strK = ToString( K );

  String in_bases = IN_HEAD + ".unibases.k" + strK;
  String in_unipaths = IN_HEAD + ".unipaths.k" + strK;
  String in_adj_graph = IN_HEAD + ".unipath_adjgraph.k" + strK;

  String out_log = OUT_HEAD + ".SplitUnibases.log";
  String out_bases = OUT_HEAD + ".fastb";
  String out_paths = OUT_HEAD + ".paths.k" + strK;
  String out_pairs = OUT_HEAD + ".pairs";

  // Load.

  cout << Date( ) << ": loading unibases" << endl;
  vecbvec bases( in_bases );
  cout << Date( ) << ": loading unipaths" << endl;
  vecKmerPath paths( in_unipaths );

  // Add junctions.

  if (USE_ADJ_GRAPH) {
    cout << Date( ) << ": loading unipath adj graph" << endl;
    digraph unigraph;
    BinaryRead( in_adj_graph, unigraph );

    cout << Date( ) << ": adding unipath links" << endl;
    int nuni = bases.size( );
    vecbasevector all(bases);
    vecKmerPath allp(paths);
    for ( int id1 = 0; id1 < nuni; id1++ ) {
      vec<int> nexts = unigraph.From(id1);
      for ( size_t j = 0; j < nexts.size( ); j++ ) {
    	int id2 = nexts[j];
    	basevector b = bases[id1];
    	b.resize( b.size( ) + 1 );
    	b.Set( b.size( ) - 1, bases[id2][K-1] );
    	all.push_back_reserve(b);
    	KmerPath p = paths[id1];
    	 p.Append( paths[id2] );
    	 allp.push_back_reserve(p);
      } 
    }
    bases = all;
    paths = allp;
    
  } else {
    cout << Date( ) << ": computing unipath adj from k-1 mers" << endl;
    vec< vec<int> > nexts;
    GetNexts( K, bases, nexts );

    cout << Date( ) << ": adding unipath links" << endl;
    int nuni = bases.size( );
    vecbasevector all(bases);
    vecKmerPath allp(paths);
    for ( int id1 = 0; id1 < nuni; id1++ )
    {    for ( int j = 0; j < nexts[id1].isize( ); j++ )
         {    int id2 = nexts[id1][j];
              basevector b = bases[id1];
	      b.resize( b.size( ) + 1 );
	      b.Set( b.size( ) - 1, bases[id2][K-1] );
	      all.push_back_reserve(b);
	      KmerPath p = paths[id1];
	      p.Append( paths[id2] );
	      allp.push_back_reserve(p);    }    }
    bases = all;
    paths = allp;

  }
  // Open streams.

  ofstream log( out_log.c_str( ) );
  IncrementalWriter<bvec> bases_out( out_bases.c_str( ) );
  IncrementalWriter<KmerPath> paths_out( out_paths.c_str( ) );

  // Parse fastb vector.

  int dotter = 100000;
  cout << Date( ) << ": parsing "
       << bases.size( ) << " bvecs (.="
       << dotter << " obejcts)"
       << endl;

  unsigned int n_reads = 0;
  for (size_t read_id=0; read_id<bases.size( ); read_id++) {
    if ( read_id % dotter == 0 ) Dot( cout, read_id / dotter );

    // This should never happen.
    const bvec &read = bases[read_id];
    const KmerPath &path = paths[read_id];
    if ( (int)read.size( ) < K ) {
      log << "WARNING: read_" << read_id
	  << " has length " << read.size( )
	  << ", which is less than K (=" << K
	  << ")\n";
      continue;
    }

    // Loop over chunks.
    unsigned int cursor = 0;
    while ( cursor < read.size( ) ) {
      n_reads++;
      int end = Min( cursor + (unsigned int)CHUNK_LEN, read.size( ) );
      bvec chunk( read, cursor, end - cursor );
      log << "r_" << read_id << " [" << cursor << ", " << end << ")\n";

      // Bases.
      bases_out.add( chunk );

      // vecKmerPath.
      KmerPathLoc start_kpl = path.BasePosToSegmentPos( cursor );
      KmerPathLoc stop_kpl = path.BasePosToSegmentPos( end - K );
      KmerPath chunk_path;

      path.CopySubpath( start_kpl, stop_kpl, chunk_path );
      paths_out.add( chunk_path );

      // Update cursor.
      if ( end == (int)read.size( ) ) break;
      cursor = end - K;
    }
  }
  cout << endl;

  // Close streams.
  cout << Date( ) << ": closing streams" << endl;
  bases_out.close( );
  paths_out.close( );
  log.close( );

  // Save empty pairs file.
  cout << Date( ) << ": saving (empty) pairs file" << endl;

  PairsManager emptypairs( n_reads );
  emptypairs.Write( out_pairs );

  // Done.
  cout << Date( ) << ": done" << endl;

}
