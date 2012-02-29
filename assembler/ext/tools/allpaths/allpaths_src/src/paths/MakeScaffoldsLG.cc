///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "CoreTools.h"
#include "MainTools.h"

#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "pairwise_aligners/AlignConsecutiveContigs.h"
#include "paths/Alignlet.h"
#include "paths/FixScaffoldsCore.h"
#include "paths/IdentifyCircularScaffolds.h"
#include "paths/MakeScaffoldsCloseBest.h"
#include "paths/MakeScaffoldsScoredBest.h"
#include "paths/RegapSupers.h"
#include "paths/SaveScaffoldGraph.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/PairDistCorrection.h"
#include "paths/PdfEntry.h"
#include "paths/reporting/ReftigUtils.h"
#include "paths/UnibaseUtils.h"

// MakeDepend: dependency QueryLookupTable
// MakeDepend: dependency MakeLookupTable


/**
 * MakeScaffoldsLG
 *
 * Scaffolding code for ALLPATHS-LG.
 * 
 * Builds scaffold graphs and merges/fixes the scaffolds using two
 * different algorithms (MakeScaffoldsScoredBest: join scaffolds by
 * using first the best edge-scores; and MakeScaffoldsCloseBest: join
 * scaffolds by looking at the most promising extension).
 *
 * HEAD_INTERIM: if not empty, save interim dir specified by STEP_INTERIM
 * STEP_INTERIM: eg "14.2", save interim step 2 of iteration MIN_LINKS=14
 * ALIGNS_OUT: if not empty, set the final set of alignments
 * MIN_LINKS: min links to use in each iterative step
 * MAX_OVERLAP: do not merge scaffolds overlapping more than this
 * MIN_REACH_AWAY: argument of FixScaffoldsCore
 * FIX_SCAFFOLDS: run FixScaffolds between merging steps
 * SUCK_SCAFFOLDS: run SuckScaffold to absorb small scaffolds in larger ones
 * REGAP: regap scaffolds
 * REGAP_NEG_GAPS: get rid of illogical negative gaps (in RegapSupers)
 * USE_HIGH_CN_ALIGNS: do a single merging run at the end with high cn aligns
 * MERGE_CONSECUTIVE_CONTIGS: try to merge consecutive, overlapping contigs
 * USE_UNIBASES: align high copy number unibases to scaffolds and break if no o
 *               no spanning read pair present
 * VERBOSE: verbose log
 */

int main( int argc, char *argv[] )
{

  RunTime( );

  BeginCommandArguments;
  CommandArgument_String( PRE );
  CommandArgument_String( DATA );
  CommandArgument_String( RUN );
  CommandArgument_String_OrDefault( SUBDIR, "test" );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_String_OrDefault( READS, "scaffold_reads" );
  CommandArgument_String_OrDefault( ALIGNS, "scaffold_reads_filtered" );
  CommandArgument_String_OrDefault( SCAFFOLDS_IN, "initial_scaffolds" );
  CommandArgument_String_OrDefault( HEAD_OUT, "final_scaffolds" );
  CommandArgument_String_OrDefault( HEAD_INTERIM, "" );
  CommandArgument_String_OrDefault( STEP_INTERIM, "" );
  CommandArgument_String_OrDefault( ALIGNS_OUT, "" );
  CommandArgument_String_OrDefault( MIN_LINKS, "{14,6,3,2}" );
  CommandArgument_Int_OrDefault( MAX_OVERLAP, 15000 );
  CommandArgument_Int_OrDefault( MIN_REACH_AWAY, 4 );
  CommandArgument_Bool_OrDefault( FIX_SCAFFOLDS, True );
  CommandArgument_Bool_OrDefault( SUCK_SCAFFOLDS, True );
  CommandArgument_Bool_OrDefault( REGAP, True );
  CommandArgument_Bool_OrDefault( REGAP_NEG_GAPS, True );
  CommandArgument_Bool_OrDefault( USE_HIGH_CN_ALIGNS, False );
  CommandArgument_Bool_OrDefault( MERGE_CONSECUTIVE_CONTIGS, False);
  CommandArgument_Bool_OrDefault( FIT_GAPS, False );

  CommandArgument_Bool_OrDefault_Doc( USE_UNIBASES, True,
    "Use linking information over repetetive unibases to fix assemblies in FixScaffolds");
  CommandArgument_String_OrDefault_Doc( UNIBASES_HEAD, "all_reads",
    "Unibases to be used in extended FixScaffolds algorithm");
  CommandArgument_Int_OrDefault_Doc( UNIBASES_K, 96,
    "K-mer size used for generating unibases (unipaths)");

  CommandArgument_Bool_OrDefault( VERBOSE, False );
  EndCommandArguments;

  // Thread control (OMP used in FixScaffoldsCore)
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Dir and file names.
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;
  String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
  String tmp_dir = sub_dir + "/tmp";

  String pairs_file = run_dir + "/" + READS + ".pairs";
  String aligns_file = sub_dir + "/" + ALIGNS + ".qltoutlet";
  String index_file =  sub_dir + "/" + ALIGNS + ".qltoutlet.index";
  String unfilt = USE_HIGH_CN_ALIGNS ? ALIGNS.Before( "_filtered" ) : "";
  String unfilt_index_file = sub_dir + "/" + unfilt + ".qltoutlet.index";
  String contigs_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
  String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
  
  String unibases_file = 
    run_dir + "/" + UNIBASES_HEAD + ".unibases.k" + ToString(UNIBASES_K);
  String unibases_cn_file = 
    run_dir + "/" + UNIBASES_HEAD + ".unipaths.predicted_count.k" + ToString(UNIBASES_K);

  String final_head = sub_dir + "/" + HEAD_OUT;
  

  // Find ploidy.
  const int ploidy = FirstLineOfFile( run_dir + "/ploidy" ).Int( );
  ForceAssertGt( ploidy, 0 );

  Mkpath( tmp_dir );
  
  // Load.
  vec<int> min_links;
  ParseIntSet( MIN_LINKS, min_links, false );

  cout << Date( ) << ": loading contigs fasta" << endl;
  vec<fastavector> contigs;
  LoadFromFastaFile( contigs_file, contigs );
  
  cout << Date( ) << ": loading supers" << endl;
  vec<superb> supers;
  ReadSuperbs( supers_file, supers );
  
  cout << Date( ) << ": loading index" << endl;
  vec<int> index;
  BinaryRead3( index_file, index );
  
  vec<int> unfilt_index;
  if ( USE_HIGH_CN_ALIGNS ) {
    cout << Date( ) << ": loading unfiltered index" << endl;
    BinaryRead3( unfilt_index_file, unfilt_index );
    ForceAssertEq(unfilt_index.size(), index.size());
  }

  cout << Date( ) << ": loading aligns" << endl;
  vec<alignlet> aligns;
  BinaryRead3( aligns_file, aligns );
  
  cout << Date( ) << ": loading pairs" << endl;
  PairsManager pairs( pairs_file );
  size_t n_pairs = pairs.nPairs( );
  
  // pair length correction 
  // correction
  vec<int> seps_empty;
  vec<int> seps;
  vec<int> sds;
  if (FIT_GAPS) {
    String reads_head =  run_dir + "/" + READS;
    PairDistCorrectionFromIntDistOld(reads_head, pairs, aligns, index, seps, sds, VERBOSE);
  }
  cout << Date( ) << ": done loading\n" << endl;

  // These are needed to run FixScaffoldsCore.
  vec<int> trace_ids;    
  ofstream devnull ( "/dev/null" );
  ostream &fs_log = VERBOSE ? cout : devnull;

  
  vec<alignlet> highCNaligns; // unibase alignment data
  vec< vec<int> > highCN_index; // unibase alignment index

  if ( USE_UNIBASES ){
    cout << Date() << ": Using unibases for scaffold fixing" << endl;
    int cn_thresh      = ploidy;  // threshold for high copy number unibase (greater than)
    prob_t prob_thresh = 0;       // copy number probability prediction threshold
    
    vecbasevector highCNunis;

    {
      vecbasevector unibases( unibases_file );
      
      vec<int> toRc;
      UnibaseInvolution( unibases, toRc, UNIBASES_K );
      vec<Bool> uused( unibases.size(), False );
      

      ForceAssert( IsRegularFile( unibases_cn_file ) );
      VecPdfEntryVec ucn_pdfs;
      ucn_pdfs.ReadAll( unibases_cn_file ); // read unibase copy number predictions
      
      ForceAssertEq( unibases.size(), ucn_pdfs.size() );
      
      for ( size_t uid = 0; uid < unibases.size(); uid++ ){
	if ( uused[uid] ) continue;
	uused[uid] = uused[ toRc[uid] ] = True;
	int cn = -1;
	prob_t maxp = 0;
	for ( unsigned int pi = 0; pi < ucn_pdfs[uid].size( ); pi++ ){    
	  if ( ucn_pdfs[uid][pi].second > maxp ) {    
	    cn = ucn_pdfs[uid][pi].first;
	    maxp = ucn_pdfs[uid][pi].second;    
	  }    
	} 
	if ( cn > cn_thresh && maxp > prob_thresh )
	  highCNunis.push_back( unibases[uid] );
      }   
    }
    cout << "Found " << highCNunis.size() << " highCN unibases" << endl;

    String work_dir = sub_dir + "/highCNunibases";
    if ( ! IsRegularFile(work_dir) )
      Mkdir777(work_dir);

    String lookup_file = work_dir + "/" + SCAFFOLDS_IN + ".contigs.lookup";
    if ( ! IsRegularFile( lookup_file ) ){
      cout << Date() << ": Creating a lookup table" << endl;
      SystemSucceed( "MakeLookupTable" + ARG(SOURCE, contigs_file )
		     + ARG(OUT_HEAD, lookup_file.SafeBefore(".lookup")) + ARG(LOOKUP_ONLY, True)
			   + ARG( NH, True ) );
    }

    cout << Date() << ": Obtaining high-CN unibase alignments" << endl;
    highCNaligns.reserve( highCNunis.size() );
    highCN_index.resize( highCNunis.size() );
    
    String highCNunis_file = work_dir + "/" + UNIBASES_HEAD + ".unibases.highCN.fastb";
    highCNunis.WriteAll( highCNunis_file );

    String highCNlook_aligns_file = work_dir + "/" + UNIBASES_HEAD + ".unibases.highCN.qltout";    
    vec<look_align> highCNlook_aligns;
    SystemSucceed( "Fasta2Fastb" + ARG(IN, contigs_file )
		   + ARG(OUT, lookup_file.SafeBefore(".lookup") + ".fastb" ) + ARG( NAMES, False ) +
			   ARG( NH, True ) );
    GetAlignsFast( 12, highCNunis_file, lookup_file, highCNlook_aligns_file, highCNlook_aligns, false, work_dir );
   
    cout << "Found " << highCNlook_aligns.size() << " unibase alignments.\nCreating index." << endl;
    for ( int ai = 0; ai < highCNlook_aligns.isize(); ai++ ){
      size_t qid = highCNlook_aligns[ai].QueryId();
      highCNaligns.push_back( highCNlook_aligns[ai] );
      highCN_index[qid].push_back( ai );
    }
    int nAlignsNotFound = 0; 
    vec<int> unalignedLens;
    for ( size_t qid = 0; qid < highCN_index.size(); qid++ ){
      if ( highCN_index[qid].size() == 0 ){
	nAlignsNotFound++;
	unalignedLens.push_back( highCNunis[qid].size() );
      }
    }
    cout << "INFO: did not find alignments for " << nAlignsNotFound << " high-CN unibases" << endl;
    cout << " their lengths are: "; unalignedLens.Print(cout); cout << endl;
    
    if ( ! VERBOSE ){
      // clean up
      Remove(lookup_file);
      Remove( lookup_file.SafeBefore(".lookup") + ".fastb" );

    }
    
  }

  // Parse STEP_INTERIM.
  int links_interim = -1;
  int iter_interim = -1;
  if ( HEAD_INTERIM != "" ) {
    links_interim = STEP_INTERIM.Before( "." ).Int( );
    iter_interim = STEP_INTERIM.After( "." ).Int( );
  }

  // Keep track if scaffolds tagged as circular.
  vec<Bool> is_circular;

  // Main loop over min_links.
  for (int jj=0; jj<min_links.isize( ); jj++) {
    if ( supers.size( ) < 2 ) break;

    int nlinks = min_links[jj];
    cout << Date( ) << ": starting loop with min_links = " << nlinks << endl;
    ReportScaffoldsBrief( supers, nlinks, 0, cout );
    if ( links_interim == nlinks && iter_interim == 0 ) {
      String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_0";
      SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	  supers, &aligns, &index );
    }

    // Adjust MAX_OVERLAP.
    double adjust_factor = double( jj + 1 ) / double( min_links.size( ) );
    int max_overlap = int ( adjust_factor * double( MAX_OVERLAP ) );

    // Scored best and fix.
    int minsep = - max_overlap;
    MakeScaffoldsScoredBest( pairs, contigs, supers, aligns, index, cout,
	&nlinks, &minsep, SUCK_SCAFFOLDS, VERBOSE );
    ReportScaffoldsBrief( supers, nlinks, 1, cout );
    if ( links_interim == nlinks && iter_interim == 1 ) {
      String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_1";
      SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	  supers, &aligns, &index );
    }

    if ( FIX_SCAFFOLDS )
      if (USE_HIGH_CN_ALIGNS)
	FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs,
			  supers, aligns, index, unfilt_index, highCNaligns, highCN_index, fs_log, VERBOSE );
      else
	FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs,
			  supers, aligns, index, index, highCNaligns, highCN_index, fs_log, VERBOSE );
    ReportScaffoldsBrief( supers, nlinks, 2, cout );
    if ( links_interim == nlinks && iter_interim == 2 ) {
      String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_2";
      SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	  supers, &aligns, &index );
    }

    IdentifyCircularScaffolds( pairs, contigs, supers, aligns, index,
	is_circular, ( VERBOSE ? &cout : 0 ), 1 );

    // Close best and fix.
    String str_links = "{" + ToString( nlinks ) + "}";
    MakeScaffoldsCloseBest( supers, contigs, aligns, index,
	pairs, "", "", str_links, "", max_overlap );
    ReportScaffoldsBrief( supers, nlinks, 3, cout );
    if ( links_interim == nlinks && iter_interim == 3 ) {
      String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_3";
      SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	  supers, &aligns, &index );
    }

    if ( FIX_SCAFFOLDS )
      if (USE_HIGH_CN_ALIGNS)
	FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs,
			  supers, aligns, index, unfilt_index, highCNaligns, highCN_index, fs_log, VERBOSE );
      else
	FixScaffoldsCore( trace_ids, pairs, MIN_REACH_AWAY, contigs,
			  supers, aligns, index, index, highCNaligns, highCN_index, fs_log, VERBOSE );
    ReportScaffoldsBrief( supers, nlinks, 4, cout );
    if ( links_interim == nlinks && iter_interim == 4 ) {
      String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_4";
      SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	  supers, &aligns, &index );
    }

    IdentifyCircularScaffolds( pairs, contigs, supers, aligns, index,
	is_circular, ( VERBOSE ? &cout : 0 ), 1 );

    // Regap.
    if ( REGAP ) {
      // Apply the gap correction if FIT_GAPS switch turned on
      if ( FIT_GAPS)
	pairs.AddSeps(seps);
      RegapSupers( fs_log, supers, pairs, aligns, index, VERBOSE, max_overlap,
	  1, 12000, REGAP_NEG_GAPS );
      if ( FIT_GAPS)
	pairs.AddSeps(seps_empty);
      ReportScaffoldsBrief( supers, nlinks, 5, cout );
      if ( links_interim == nlinks && iter_interim == 5 ) {
	String head_out = HEAD_INTERIM + ".ml_" + ToString( nlinks ) + "_5";
	SaveInterimScaffolds( data_dir, head_out, pairs, contigs, 
	    supers, &aligns, &index );
      }
    }

    cout << endl;

  }

  // Run one more merging iteration (optional).
  if ( USE_HIGH_CN_ALIGNS ) {
    cout << Date( ) << ": starting extra loop (high CN aligns)" << endl;

    int min_links = 2;
    int max_links = 28;
    int min_scaffold_len = 100000;
    ReportScaffoldsBrief( supers, min_links, 0, cout );

    MakeScaffoldsCloseBest( supers, contigs, aligns, unfilt_index, pairs, "",
	"", ToString( min_links ), ToString( max_links ),
	MAX_OVERLAP, min_scaffold_len );
    ReportScaffoldsBrief( supers, min_links, 2, cout );
    cout << endl;
  }

  // Merge consecutive contigs.
  if ( MERGE_CONSECUTIVE_CONTIGS ) {
    const int MIN_OVERLAP = 0;
    const float DEV_FLEX = 3.0;

    vecbvec bcontigs;
    bcontigs.reserve( contigs.size( ) );
    for (int ii=0; ii<contigs.isize( ); ii++)
      bcontigs.push_back( contigs[ii].ToBasevector( ) );

    int dotter = 100;
    cout << Date( ) << ": merging consecutive contigs (.="
      << dotter << " supers, with "
      << supers.isize( ) << " supers in input):" << endl;

    for (int super_id=0; super_id<supers.isize( ); super_id++) {
      if ( super_id % dotter == 0 ) Dot( cout, super_id / dotter );
      const superb &sup = supers[super_id];
      for (int cgpos=0; cgpos<sup.Ntigs( )-1; cgpos++) {
	int gap = sup.Gap( cgpos );
	int dev = sup.Dev( cgpos );
	int flex_gap = gap - int( DEV_FLEX * double( dev ) );
	if ( flex_gap <= - MIN_OVERLAP ) {
	  AlignConsecutiveContigs( super_id, cgpos, bcontigs, supers );
	}
      }
    }
    cout << endl;

    vec<fastavector> new_contigs;
    new_contigs.reserve( bcontigs.size( ) );
    for (size_t ii=0; ii<bcontigs.size( ); ii++)
      new_contigs.push_back( fastavector( bcontigs[ii] ) );

    swap( contigs, new_contigs );
  }

  // Save ids of scaffolds tagged by IdentifyCircularScaffolds.
  cout << Date( ) << ": saving circular scaffolds info" << endl;
  String out_file = final_head + ".is_circular";
  IdentifyCircularScaffolds( pairs, contigs, supers, aligns, index,
      is_circular, ( VERBOSE ? &cout : 0 ), 1 );
  BinaryWrite3( out_file, is_circular );

  // Save aligns and index.
  if ( ALIGNS_OUT != "" ) {
    cout << Date( ) << ": saving aligns and index" << endl;
    String out_aligns_file = sub_dir + "/" + ALIGNS_OUT + ".qltoutlet";
    String out_index_file =  sub_dir + "/" + ALIGNS_OUT + ".qltoutlet.index";
    BinaryWrite3( out_aligns_file, aligns );
    BinaryWrite3( out_index_file, index );
  }

  // Save and leave.
  cout << Date( ) << ": saving final scaffolds" << endl;
  SaveScaffoldAssembly( final_head, supers, contigs );
  WriteSummary( final_head + ".summary", supers);

  cout << Date( ) << ": done" << endl;

}

