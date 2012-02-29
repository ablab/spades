/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2010) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include <omp.h>

#include "MainTools.h"
#include "Basevector.h"
#include "Fastavector.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "PairsManager.h"
#include "VecTemplate.h"
#include "lookup/FirstLookupFinder.h"
#include "lookup/LookAlign.h"
#include "lookup/LookupTabBuilder.h"
#include "paths/Alignlet.h"
#include "paths/CLinFastavec.h"
#include "paths/RemoveDuplicateAligns.h"
#include "paths/SquashHyperFastavector.h"
#include "util/SearchFastb2Core.h"

/**
 * AlignPairsToFasta
 *
 * AlignsPairs to a fasta file representing the target (e.g.,
 * contigs).  Uses a perfect aligner.  Saves the alignments and a
 * per-read index. --bruce 3/31/10
 *
 * ALIGNS: it defaults to READS
 * REMOVE_DUPLICATE: if true, run RemoveDuplicateAligns
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
  CommandArgument_String( FASTA );
  CommandArgument_Int_OrDefault_Doc( K_MIN, 12, "minimum k-mer size to use in lookup" );
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "" );
  CommandArgument_String_OrDefault( WRUN, "recover" );
  CommandArgument_Bool_OrDefault( REMOVE_DUPLICATES, True );
  CommandArgument_Bool_OrDefault( VERBOSE, True );
  EndCommandArguments;

  // Default args.
  if ( ALIGNS == "" ) ALIGNS = READS;

  // Thread control (Uses OMP in SearchFastb)
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String run_dir = PRE + "/" + DATA + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String wrun_dir = sub_dir + "/" + WRUN;
  
  String reads_file = run_dir + "/" + READS + ".fastb";
  String pairs_file = run_dir + "/" + READS + ".pairs";
  String target_file = sub_dir + "/" + FASTA + ".fasta";

  String out_qltlet_file = sub_dir + "/" + ALIGNS + ".qltoutlet"; 
  String out_qltlet_index_file = out_qltlet_file + ".index"; 

  // Load reads and kill the unpaired ones.

  cout << Date( ) << ": loading reads" << endl;
  vecbvec reads( reads_file );
  cout << Date( ) << ": loading pairing info" << endl;
  PairsManager pairs( pairs_file );
  {
    for ( size_t p = 0; p < reads.size( ); p++ )
      if ( pairs.isUnpaired(p) ) reads[p].clear().shrink_to_fit();
  }
  
  cout << Date( ) << ": loading target fasta" << endl;
  vecbasevector target;
  FetchReads( target, 0, target_file );

  cout << Date( ) << ": converting target to fastb" << endl;
  Mkdir777( sub_dir + "/tmp" );
  String target_fastb = sub_dir + "/tmp/AlignPairsToFasta.target.fastb";
  target.WriteAll( target_fastb );


  // Align reads.
  cout << Date( ) << ": aligning reads to fasta" << endl;
  
  vec<int> allowedKs("{12,20,26,40,60,64,80,88,96,128,144}"); // values from SearchFastb2Core.cc
  Sort( allowedKs );
  //cout << "allowedKs "; allowedKs.Print(cout); cout << endl;
  Bool foundK = False;
  for ( size_t i = 0; i < allowedKs.size(); i++ )
    if ( allowedKs[i] >= K_MIN ){
      foundK = True;
      allowedKs.SetToSubOf( allowedKs, i, allowedKs.size() - i );
      break;
    }
  //cout << "allowedKs "; allowedKs.Print(cout); cout << endl;
  if ( ! foundK ){
    cout << "SearchFastb2: Not implemented for K_MIN=" << K_MIN << "." << endl;
    cout << "Abort." << endl;
    exit(1);
  }
  int minPosReadLen = 0;
  for ( size_t id = 0; id < reads.size(); id++ )
    if ( reads[id].isize() > 0 ){
      minPosReadLen = reads[id].isize();
      break;
    }
  for ( size_t id = 0; id < reads.size(); id++ )
    if ( reads[id].isize() < minPosReadLen && reads[id].size() > 0u )
      minPosReadLen = reads[id].isize();
  //  then find largest allowed k-mer 
  int Kl = allowedKs.front(); // they were sorted first
  //cout << "initial Kl =" << Kl << endl;
  for ( size_t i = 1; i < allowedKs.size(); i++ )
    if ( allowedKs[i] <= minPosReadLen )
      Kl = allowedKs[i];
  //PRINT2( minPosReadLen, Kl );
  cout << "search k-mer size = " << Kl << endl;

  String reads_file0 = sub_dir + "/tmp/AlignPairsToFasta.reads.fastb";
  reads.WriteAll(reads_file0);
  vec< triple<int64_t,int64_t,int> > aligns;
  const int max_placements = 1;
  SearchFastb2( reads_file0, target_fastb, Kl, &aligns, 0, max_placements, 0.90, VERBOSE );
  Remove(reads_file0);
  Remove(target_fastb);

  vec<alignlet> aligns0;
  aligns0.reserve( aligns.size( ) );
  for ( size_t i = 0; i < aligns.size( ); i++ ) {
    int id = aligns[i].first, t = aligns[i].second;
    int pos2 = ( aligns[i].third >= 0 ? aligns[i].third : -aligns[i].third - 1 );
    aligns0.push( pos2, pos2 + reads[id].isize( ),
		  t, target[t].size( ), aligns[i].third >= 0 );
  }

  cout << Date( ) << ": generating aligns index" << endl;
  vec<int> aligns0_index( reads.size( ), -2 );
  for ( size_t i = 0; i < aligns.size( ); i++ ) {
    int id = aligns[i].first;
    if ( aligns0_index[id] == -2 ) aligns0_index[id] = i;
    else aligns0_index[id] = -1;
  }

  // Remove duplicates.
  
  if ( REMOVE_DUPLICATES )
    RemoveDuplicateAligns( pairs, aligns0, aligns0_index, cout );

  // Save aligns.

  cout << Date( ) << ": saving " << aligns0.size( ) << " aligns" << endl;
  BinaryWrite3( out_qltlet_file, aligns0 );
  BinaryWrite3( out_qltlet_index_file, aligns0_index );

  // Done.
  cout << Date( ) << ": done" << endl;
  _exit(0);

}
