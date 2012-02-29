///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "paths/HyperBasevector.h"
#include "paths/HyperFastavector.h"
#include "paths/FlattenHyperFastavector.h"

/**
 * FlattenHKP
 *
 * Convert input HyperKmerPath into a vec<fastavector>.  (ie flatten
 * it: this is done to improve sensitivity in the alignment step).
 * Convert to scaffold structure.
 * --bruce 3/31/10
 *
 * HYPER: input HyperKmerPath, relative to SUBDIR
 * SCAFFOLDS: output (.superb, .contigs.fasta, .assembly.fasta)
 * WRUN: unipaths info (relative to SUBDIR: needed by the KmerBaseBroker)
 * SAVE_INTERMEDIATE: save intermediete (iterative) steps from linearization
 * NUM_THREADS: number of threads for FirstLookupFinder
 */

int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String( SUBDIR );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( HYPER, "hyper_plus" );
  CommandArgument_String_OrDefault( SCAFFOLDS, "initial_scaffolds" );
  CommandArgument_String_OrDefault( WRUN, "recover" );
  CommandArgument_Int_OrDefault( MAX_CELL_SIZE, 20 );
  CommandArgument_Int_OrDefault( MIN_EDGE_TO_SAVE, 1000 );
  CommandArgument_Bool_OrDefault( INITIAL_SCAFFOLD_PER_CONTIG, True );
  CommandArgument_Bool_OrDefault( SAVE_INTERMEDIATE, False );
  CommandArgument_Bool_OrDefault( DUMP_HFV, False );
  CommandArgument_Bool_OrDefault( WRITE, True );
  CommandArgument_Bool_OrDefault( NEW_ALGORITHM, False );
  EndCommandArguments;

  // Thread control 

  NUM_THREADS = configNumThreads(NUM_THREADS);
    
  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String wrun_dir = sub_dir + "/" + WRUN;
  
  String hyper_file = sub_dir + "/" + HYPER;
  String base_out = WRITE ? sub_dir + "/" + SCAFFOLDS : "";
  String hyper_inter_file = SAVE_INTERMEDIATE ? hyper_file + "_inter" : "";
  String dump_hfv_head  = DUMP_HFV ? hyper_file + ".fastavector" : "";
  String log_file = hyper_file + ".FlattenHKP.log"; 
  
  HyperFastavector hfv;

  // The log stream. 
  cout << "\nSending log to: " << log_file << "\n" << endl;
  ofstream log( log_file.c_str( ) );
  PrintCommandPretty( log );

  log << Date( ) << ": Loading HyperKmerPath" << endl;
  HyperKmerPath hk( hyper_file );
  int K = hk.K( );

  log << Date( ) << ": Loading KmerBaseBroker" << endl;
  KmerBaseBroker kbb( wrun_dir, K );

  log << Date( ) << ": Turning HyperKmerPath into HyperFastavector" << endl;
  hfv = HyperFastavector( HyperBasevector( hk, kbb ) );
  
  // Flatten hyper.

  FlattenHyperFastavector( log,
			   hfv,
			   base_out,
			   hyper_inter_file,
			   dump_hfv_head,
			   NEW_ALGORITHM,
			   INITIAL_SCAFFOLD_PER_CONTIG,
			   MAX_CELL_SIZE,
			   MIN_EDGE_TO_SAVE,
			   NUM_THREADS );
  
  // Done.
  
  log << Date( ) << ": done" << endl;
  return 0;
  
}
