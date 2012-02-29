///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include <algorithm>
#include "CoreTools.h"
#include "MainTools.h"
#include "Fastavector.h"
#include "Superb.h"
#include "TokenizeString.h"
#include "paths/SaveScaffoldGraph.h"

/**
 * AgpToAllPaths
 *
 * Load agp and fasta file and generate an AllPaths assembly.
 *
 * HEAD_IN: it loads <HEAD_IN>.{fasta,agp}
 * HEAD_OUT: output is saved with SaveScaffoldAssembly
 * DEV_MULTIPLIER: set dev to ( DEV_MULTIPLIER * gap ), where 0 is allowed
 */
int main( int argc, char *argv[] )
{ 
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( HEAD_IN );
  CommandArgument_String( HEAD_OUT );
  CommandArgument_Double_OrDefault( DEV_MULTIPLIER, 0.25 );
  CommandArgument_Bool_OrDefault( VERBOSE, False );
  CommandArgument_Bool_OrDefault( SAVE_FASTB, False );
  EndCommandArguments;

  // File names.
  String agp_file = HEAD_IN + ".agp";
  String fasta_file = HEAD_IN + ".fasta";

  // Check for files.
  if ( ! IsRegularFile( agp_file ) ) {
    cout << "Fatal error, agp file not found.\nExit.\n" << endl;
    return 1;
  }
  if ( ! IsRegularFile( fasta_file ) ) {
    cout << "Fatal error, fasta file not found.\nExit.\n" << endl;
    return 1;
  }

  // Load.
  cout << Date( ) << ": loading agp" << endl;
  vec< vec<String> > agp;
  String line;
  vec<String> tokens;
  ifstream in( agp_file.c_str( ) );
  while ( in ) {
    getline( in, line );
    if ( ! in ) break;
    Tokenize( line, tokens );
    agp.push_back( tokens );
  }
  
  cout << Date( ) << ": loading fasta" << endl;
  vec<fastavector> contigs;
  vec<String> vnames;
  LoadFromFastaFile( fasta_file, contigs, vnames );

  // Generate superb.
  cout << Date( ) << ": generating super structure" << endl;
  String previous_super_name = "";
  vec<superb> supers;
  for (size_t ii=0; ii<agp.size( ); ii++) {
    bool is_gap = agp[ii][6].Contains( "fragment", 0 );
    if ( is_gap ) {
      int gap_size = agp[ii][5].Int( );
      int gap_dev = Max( 0, int( double( gap_size ) * DEV_MULTIPLIER ) );
      supers[ supers.size( )-1 ].AppendGap( gap_size, gap_dev );
      continue;
    }

    size_t contig_id = find(vnames.begin(), vnames.end(), agp[ii][5] ) - vnames.begin();
    if (contig_id == vnames.size()) {
      cout << "Missing contig in fasta file: " << agp[ii][5] << endl;
      FatalErr("Mismatch between fasta and AGP files");
    }

    int contig_len = agp[ii][7].Int( );

    String super_name = agp[ii][0];
    if (super_name != previous_super_name) {
      superb sup;
      sup.PlaceFirstTig( contig_id, contig_len );
      supers.push_back( sup );
      previous_super_name = super_name;
    }
    else
      supers[ supers.size( )-1 ].AppendTig( contig_id, contig_len );
  }

  // Save.
  cout << Date( ) << ": saving assembly" << endl;
  ostream *log = VERBOSE ? &cout : 0;
  SaveScaffoldAssembly( HEAD_OUT, supers, contigs, log, SAVE_FASTB );
  WriteSummary( HEAD_OUT + ".summary", supers );


  // Done.
  cout << Date( ) << ": done" << endl;
  
}
