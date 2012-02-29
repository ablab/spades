///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "MainTools.h"

#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/IdentifyCircularScaffolds.h"

/**
 * TagCircularScaffolds
 *
 * Run IdentifyCircularScaffolds to identify and tag circular
 * scaffolds.
 *
 * NUM_THREADS: if negative, use all available
 */
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_String_OrDefault( FRAG_READS, "" );
  CommandArgument_String_OrDefault( JUMP_READS, "" );
  CommandArgument_String_OrDefault( LONG_JUMP_READS, "" );
  CommandArgument_String_OrDefault( ASSEMBLY, "" );
  CommandArgument_UnsignedInt_OrDefault_Doc( NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;

  ForceAssert(FRAG_READS != "" | JUMP_READS != "" | LONG_JUMP_READS != "");

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

  String head = sub_dir + "/" + ASSEMBLY;
  String contigs_file = head + ".contigs.fasta";
  String supers_file = head + ".superb";
  String rings_file = head + ".rings";
  
  // Load.
  read_locs_on_disk locs_file( head, run_dir );

  cout << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );

  // Import pairs
  vec<String> pairs_files(3, "");
  vec<PairsManager> pairs;
  if (FRAG_READS != "")
    pairs_files[0] = FRAG_READS;
  if (JUMP_READS != "")
    pairs_files[1] = JUMP_READS;
  if (LONG_JUMP_READS != "")
    pairs_files[2] = LONG_JUMP_READS;
 
  for (int ii = 0; ii < 3; ii++) {
    if ( pairs_files[ii] != "" ) {
      if ( IsRegularFile( run_dir + "/" + pairs_files[ii] + ".pairs" ) ) {
	cout << Date( ) << ": loading " << pairs_files[ii] << " pairs file" << endl;
	pairs.push_back( PairsManager( run_dir + "/" + pairs_files[ii] + ".pairs" ) );
      } else {
	FatalErr("Unable to find pairs file for: " + pairs_files[ii]);
      }
    } else
      pairs.push_back( PairsManager( ) );
  }
  
  cout << Date( ) << ": done loading\n" << endl;
  
  // Run code.
  vec<CRing> rings;
  IdentifyCircularScaffolds( pairs, contigs, supers, locs_file, rings,
			     NUM_THREADS, &cout );

  // Save output.
  cout << "\n" << Date( ) << " saving rings file" << endl;
  ofstream out( rings_file.c_str( ) );
  if ( rings.size( ) < 1 )
    out << "No circular scaffolds found.\n" << endl;
  else {
    vec< vec<String> > table;
    vec<String> line;
    rings[0].GenerateLegend( line );
    table.push_back( line );
    for (int ii=0; ii<rings.isize( ); ii++) {
      rings[ii].PrintableInfo( line );
      table.push_back( line );
    }
    PrintTabular( out, table, 3, "rrrr" );
    out.close( );
  }
  
  // Done.
  cout << Date( ) << ": done" << endl;
  
}

