///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include "MainTools.h"
#include "FetchReads.h"
#include "FetchReadsAmb.h"
#include "ParseSet.h"
#include "SupersHandler.h"
#include "lookup/LookAlign.h"
#include "paths/reporting/CGapStats.h"
#include "paths/reporting/EvalScaffoldsUtils.h"
#include "paths/reporting/ReftigUtils.h"
#include "util/RunCommand.h"
#include "reporting/PerfStat.h"

/**
 * EvalScaffolds
 *
 * Align <SCAFFOLDS><CONTIGS>.fasta onto a reference, and evaluate the
 * result. Output is saved in a new directory <SCAFFOLDS>.Eval
 *
 * LOOKUP: <LOOKUP>.Before( ".lookup" ) + ".fastb" is needed for genome size
 * SCAFFOLDS: it loads <SCAFFOLDS>{<CONTIGS>.fasta,.superb}
 * BRIEF: print summary only
 * FORCE: do not use cached aligns
 */
int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( LOOKUP );
  CommandArgument_String( SCAFFOLDS );
  CommandArgument_String_OrDefault( CONTIGS, ".contigs" );
  CommandArgument_String_OrDefault( SUPER_IDS, "" );
  CommandArgument_String_OrDefault( OUT_DIR, "" );
  CommandArgument_Bool_OrDefault( BRIEF, False );
  CommandArgument_Bool_OrDefault( FORCE, False );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;
  
  // Dir and file names.
  String out_dir = (OUT_DIR == "" ? SCAFFOLDS + ".Eval" : OUT_DIR);

  String genome_fastb = LOOKUP.Before( ".lookup" ) + ".fastb";

  String contigs_fasta = SCAFFOLDS + CONTIGS + ".fasta";
  String contigs_fastb = SCAFFOLDS + CONTIGS + ".fastb";
  String contigs_fastamb = SCAFFOLDS + CONTIGS + ".fastamb";
  String supers_file = SCAFFOLDS + ".superb";

  String select_file = out_dir + "/ids.select";
  String qlt_file = out_dir + "/aligns.qlt";
  String out_file = out_dir + "/main.log";

  Mkpath( out_dir );
 
  // Thread control
  // This module uses OMP indirectly, so we must set the omp thread count here.
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );
 
  // Load.
  cout << Date( ) << ": loading supers" << endl;
  shandler supers( -1, supers_file );
  
  longlong genome_size = 0;
  {
    vecbvec genome( genome_fastb );
    for (size_t ii=0; ii<genome.size( ); ii++)
      genome_size += genome[ii].size( );
  }

  // Select supers to align.
  int n_selected = 0;
  vec<bool> selected;
  if ( SUPER_IDS != "" ) {
    selected.resize( supers.Size( ), false );
    vec<int> sel;
    ParseIntSet( SUPER_IDS, sel );
    n_selected = sel.size( );
    for (int ii=0; ii<sel.isize( ); ii++)
      selected[ sel[ii] ] = true;
  }
  else {
    selected.resize( supers.Size( ), true );
  }
  
  // Extract selected contigs.
  vec<int> c_ids;
  for (int ii=0; ii<selected.isize( ); ii++) {
    if ( ! selected[ii] ) continue;
    const superb &sup = supers[ii];
    for (int pos=0; pos<sup.Ntigs( ); pos++)
      c_ids.push_back( sup.Tig( pos ) );
  }
  sort( c_ids.begin( ), c_ids.end( ) );

  if ( FORCE || ! IsRegularFile( select_file ) ) {
    cout << Date( ) << ": saving select file" << endl;
    ofstream out( select_file.c_str( ) );
    for (int ii=0; ii<c_ids.isize( ); ii++)
      out << c_ids[ii] << "\n";
    out.close( );
  }  

  // Generate contigs fastb, fastamb.
  if ( FORCE || ! IsRegularFile( contigs_fastb ) ) {
    cout << Date( ) << ": saving contigs fastb" << endl;
    vecbasevector cbases;
    FetchReads( cbases, 0, contigs_fasta );
    cbases.WriteAll( contigs_fastb );
    
    cout << Date( ) << ": saving contigs fastamb" << endl;
    vecbitvector amb;
    FetchReadsAmb( amb, contigs_fasta );
    amb.WriteAll( contigs_fastamb );
  }

  // Align selected edges, and (re-)load aligns.
  vec<look_align_plus> hits;
  if ( FORCE || ! IsRegularFile( qlt_file ) ) {
    cout << Date( ) << ": aligning " << c_ids.size( ) << " contigs" << endl;
    vec<look_align> aligns;
    // Note that this function uses OMP (badly)
    GetAlignsFast( 96, contigs_fastb, LOOKUP, qlt_file, aligns, !FORCE, out_dir ); 
  }
  
  cout << Date( ) << ": loading aligns from file\n" << endl;
  LoadLookAlignPlus( qlt_file, hits );
  
  cout << "Sending output to " << out_file << "\n" << endl;
  ofstream out( out_file.c_str( ) );
  PrintCommandPretty( out );
  
  // A map edge_id to hits.
  vec< vec<int> > edge2hits( supers.NContigs( ) );
  for (int ii=0; ii<hits.isize( ); ii++)
    edge2hits[ hits[ii].query_id ].push_back( ii );

  // Gap stats.
  CGapStats gap_stats( &supers );
  
  // Loop over all supers.
  int n_local_inversions = 0;
  vec<SPlacement> all_placs;
  for (int super_id=0; super_id<supers.Size( ); super_id++) {
    vec<int> chain;
    SuperChain( super_id, supers, hits, edge2hits, chain );
    ChainGaps( supers, hits, chain, gap_stats );
    
    vec<SPlacement> placs;
    DigestChain( supers, hits, chain, placs );
    copy( placs.begin( ), placs.end( ), back_inserter( all_placs ) );

    int n_inversions = 0;
    ofstream devnull ( "/dev/null" );
    ostream &brief_out = BRIEF ? devnull : out;
    PrintChain( supers, hits, chain, n_inversions, brief_out );

    n_local_inversions += n_inversions;
  }
  
  // Sort and print placements.
  sort( all_placs.begin( ), all_placs.end( ) );

  vec< vec<String> > table;
  for (int ii=0; ii<all_placs.isize( ); ii++) {
    vec<String> line;
    all_placs[ii].PrintInfo( supers, line );
    table.push_back( line );
  }

  if ( all_placs.size( ) > 0 ) {
    const String separator = "   ";
    all_placs[0].PrintLegend( out );
    BeautifyAndPrintTable( table, out, separator );
    out << endl;
  }
  else {
    out << "No aligns found\n" << endl;
  }
  
  // Report on gap accuracy.
  gap_stats.PrintReport( out );

  // Report on number of local misorderings.
  double genome_size_mb = double( genome_size ) / 1000000.0;
  double ratio = double( n_local_inversions ) / genome_size_mb;
  
  out << "LOCAL MISORDERINGS: " << n_local_inversions
      << " found, corresponding to " << ToString( ratio, 2 )
      << " inversions per Mb\n"
      << "(based on a genome size of " << genome_size
      << " Mb, as found in the genome.fastb file).\n" << endl;

  PerfStat::log( ) << std::fixed << std::setprecision(2)
    << PerfStat( "misordering_rate", "number of local misorderings per Mb", ratio );

  // Done.
  out << Date( ) << ": done" << endl;
  cout << Date( ) << ": done" << endl;
  
}

