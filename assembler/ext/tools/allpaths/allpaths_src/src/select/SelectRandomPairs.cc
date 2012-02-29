///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"
#include "Basevector.h"
#include "Qualvector.h"
#include "random/Shuffle.h"
#include "PairsManager.h"
#include "feudal/IncrementalWriter.h"
#include "feudal/BinaryStream.h"
#include <map>
#include "lookup/LookAlign.h"

/**
 * SelectRandomPairs
 *
 * Randomly select pairs of reads, up to a total selected read length
 * of <MB_TOTAL> Mb or fraction FRAC. A correspondence map is also
 * saved, to trace reads back to their original ids.
 *
 * READS_IN: it loads <READS_IN>.{fastb,qualb,pairs,qltout}
 * READS_OUT: it saves <READS_OUT>.{fastb,qualb,pairs,select,qltout}
 * MB_TOTAL: required total length of selected reads (Mb)
 * FRAC: fraction of total length to use (use 1.0 to select all reads)
 * SEED: seed for the random number generator
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( READS_IN );
  CommandArgument_String( READS_OUT );
  CommandArgument_UnsignedInt_OrDefault( N_PAIRS, 0 );
  CommandArgument_Double_OrDefault( MB_TOTAL, -1 );
  CommandArgument_Double_OrDefault( FRAC, -1 );
  CommandArgument_UnsignedInt_OrDefault( SEED, 666 );
  EndCommandArguments;

  // Check arguments.
  {
    ForceAssert((  (MB_TOTAL > 0) && !(FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) &&  (FRAC > 0) && !(N_PAIRS > 0) ) ||
                ( !(MB_TOTAL > 0) && !(FRAC > 0) &&  (N_PAIRS > 0) ));
  }
  
  // File names.
  String in_bases_file = READS_IN + ".fastb";
  String in_quals_file = READS_IN + ".qualb";
  String in_pairs_file = READS_IN + ".pairs";
  String in_qlt_file   = READS_IN + ".qltout";

  String out_bases_file  = READS_OUT + ".fastb";
  String out_quals_file  = READS_OUT + ".qualb";
  String out_select_file = READS_OUT + ".select";
  String out_pairs_file  = READS_OUT + ".pairs";
  String out_qlt_file    = READS_OUT + ".qltout";

  // Output dir.
  String out_dir = out_bases_file;
  if ( ! out_dir.Contains( "/" ) ) out_dir = ".";
  else out_dir = out_dir.RevBefore( "/" );
  if ( out_dir != "." ) Mkpath( out_dir );
  
  size_t n_reads = MastervecFileObjectCount( in_bases_file );
 
  // Selected read ids.
  vec<size_t> select;
  vec<longlong> maps_to( n_reads, -1 );
  
  // Pairs and bases - scoped for memory.
  {
    cout << Date( ) << ": loading pairing info" << endl;
    PairsManager pairs( in_pairs_file );
    size_t n_pairs = pairs.nPairs();
    
    // Load bases and pairs.
    cout << Date( ) << ": loading bases" << endl;
    vecbvec bases( in_bases_file );
    
    
    // Shuffle.
    cout << Date( ) << ": shuffling pairs ids" << endl;
    vec<uint64_t> shuffled;
    Shuffle64( (uint64_t)n_pairs, shuffled, (uint64_t)SEED );
  
    size_t n_keepers = 0;
    if (FRAC == 1.0 || N_PAIRS > 0) {
      n_keepers = (FRAC == 1.0 ? n_pairs : N_PAIRS);
      
      size_t tot_length = 0;
      for (size_t ii = 0; ii < n_keepers; ii++) 
        tot_length += (bases[pairs.ID1(shuffled[ii])].size() +
                       bases[pairs.ID2(shuffled[ii])].size());

      cout << endl
           << "total length of selected reads: " << tot_length << endl
           << "pairs in input:                 " << n_pairs << endl
           << "pairs selected:                 " << n_keepers << endl
           << endl;
    } 
    else {
      // Determine total Mb to keep.
      size_t required_length = 0;
      if (MB_TOTAL > 0) required_length = size_t(MB_TOTAL * 1000000.0);
      if (FRAC     > 0) required_length = uint64_t(round(double(bases.sumSizes()) * FRAC));
      
      // Determine how many pairs we want to keep.
      size_t tot_length = 0;
      for (size_t ii=0; ii<n_pairs; ii++) {
        if ( tot_length >= required_length ) break;
        tot_length += bases[ pairs.ID1( shuffled[ii] ) ].size( );
        tot_length += bases[ pairs.ID2( shuffled[ii] ) ].size( );
        n_keepers = ii+1;
      }
      cout << endl
           << "total length of selected reads: " << tot_length << endl
           << "required total length:          " << required_length << endl
           << "pairs in input:                 " << n_pairs << endl
           << "pairs selected:                 " << n_keepers << endl
           << endl;
    }
    
    // Save pairs, and populate select list.
    cout << Date( ) << ": saving pairs" << endl;

    PairsManager sel_pairs( n_keepers * 2 );
    for (size_t ii=0; ii<n_keepers; ii++) {

      size_t id1 = select.isize();
      select.push_back( pairs.ID1( shuffled[ii] ) );
      maps_to[ pairs.ID1( shuffled[ii] ) ] = id1;
      
      size_t id2 = select.isize( );
      select.push_back( pairs.ID2( shuffled[ii] ) );
      maps_to[ pairs.ID2( shuffled[ii] ) ] = id2;

      sel_pairs.addPair( id1, id2, pairs.sep( shuffled[ii] ), pairs.sd( shuffled[ii] ),
			 pairs.libraryName( shuffled[ii] ) );
    }
    sel_pairs.Write( out_pairs_file );

    cout << Date( ) << ": saving correspondence map (select)" << endl;

    BinaryWriter::writeFile( out_select_file.c_str(), select );

    cout << Date( ) << ": saving bases" << endl;
    IncrementalWriter<bvec> sel_bases(out_bases_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_bases.add( bases[select[ii]] );
    sel_bases.close();

    

  }

  // Select corresponding quals.
  if (IsRegularFile(in_quals_file)) {
    cout << Date( ) << ": loading quals" << endl;
    vecqvec quals(in_quals_file);

    cout << Date( ) << ": saving quals" << endl;
    IncrementalWriter<qvec> sel_quals(out_quals_file.c_str());
    for (size_t ii=0; ii < select.size(); ii++)
      sel_quals.add( quals[select[ii]] );
    sel_quals.close();
  }

  // Select aligns (this could be improved: no need to load all aligns).
  if ( IsRegularFile(in_qlt_file) ){
    cout << Date() << ": loading alignments" << endl;
    vec<look_align> aligns_in;
    LoadLookAligns( in_qlt_file, aligns_in );

    cout << Date( ) << ": selecting and saving alignments" << endl;
    vec<look_align> aligns_out;
    aligns_out.reserve( aligns_in.size( ) );
    for( size_t ii=0; ii<aligns_in.size( ); ii++) {
      if ( maps_to[ aligns_in[ii].query_id ] < 0 ) continue;
      look_align al = aligns_in[ii];
      al.query_id = maps_to[ aligns_in[ii].query_id ];
      aligns_out.push_back( al );
    }
    sort( aligns_out.begin( ), aligns_out.end( ) );
    
    WriteLookAligns( out_qlt_file, aligns_out );
  }

  // Done.
  cout << Date( ) << ": done" << endl;
  
}
