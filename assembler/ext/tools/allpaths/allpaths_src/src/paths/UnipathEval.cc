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

#include <list>
#include "Basevector.h"
#include "BinsVec.h"
#include "Bitvector.h"
#include "Intvector.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "PredictionStats.h"
#include "VecUtilities.h"
#include "lookup/LookAlign.h"
#include "lookup/PerfectLookup.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/reporting/ReftigUtils.h" // GetAligns
#include "paths/simulation/GenomePlacement.h"
#include "graph/Digraph.h"
#include "paths/Unipath.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnipathFixerTools.h"
#include "reporting/PerfStat.h"
#include "system/Utils.h"

typedef double gc_t;
typedef int unipath_size_t;
typedef BinsVec2 < unipath_size_t, gc_t, PredictionStats > bin2PredStat;

void EvalUnipaths ( const vecKmerPath& unipaths,
		    const vec<tagged_rpint>& unipathsdb,
		    const vecKmerPath& genome_unipaths,
		    const vec<big_tagged_rpint>& genome_unipathsdb,
		    const digraph& unipathAdjGraph,
		    ostream &out );

void EvalGraph ( const int K,
		 const String &run_dir,
		 const String &READS,
		 const String &UNIPATHS,
		 const vecKmerPath &unipaths,
		 ostream &out );

void EvalCopyNumber( const int PLOIDY,
		     const int MATRIX_LIMIT,
		     const vecbvec &unibases,
		     const VecPdfEntryVec &estimatedCopyNumber,
		     const vec<copy_num_t> &trueCopyNumber,
		     const vec<double> &gcBins,
		     const vec<int> &lenBins,
		     const vec<int> &sizes,
		     bin2PredStat &copyNumberPloidyStats,
		     ostream &out );

void AlignUnibases( const int &K,
		    const int &BADS,
		    const int &PLOIDY,
		    const int &REF_PLOIDY,
		    const int &MATRIX_LIMIT,
		    const String &GENOME,
		    const bool &SAVE_ALIGNS,
		    const bool &SHOW_MAP,
		    const bool &COPYNUM,
		    const String &data_dir,
		    const String &run_dir,
		    const vec<double> &gcBins,
		    const vec<int> &lenBins,
		    const vec<int> &sizes,
		    const vec<nkmers_t> &usizes,
		    const VecPdfEntryVec &estimatedCopyNumber,
		    temp_file &uba,
		    bin2PredStat &copyNumberPloidyStats,
		    vecbvec &unibases,
		    vec<look_align> &aligns1,
		    vec<look_align> &aligns2,
		    vec<look_align> &aligns3,
		    vec<look_align> &to_save,
		    vec<look_align> &best_aligns,
		    vec<genome_placement> &fw_placements,
		    vec<genome_placement> &rc_placements,
		    vec< pair<nkmers_t,int> > &len_id,
		    vecbitvector &cov,
		    vecbitvector &cov100,
		    vec<Bool> &perfect,
		    vec<int> &mismatches,
		    vec<int> &indels,
		    int &nseq,
		    ostream &out );



/**
   Program: UnipathEval
   
   Print statistics on unipath sizes, or on the sizes of other objects
   that look reasonably like unipaths.  If there is a reference,
   evaluate unipaths for accuracy and completeness.

   Program parameters:

     PRE - the <WGA data directory>
     DATA - the <project directory>
     RUN - the <run directory>
     
     K - the size of kmers that define the unipaths
     
     UNIPATHS - part of the file name in the <run dir> from which the
        unipaths are read; the full name will be reads.UNIPATHS.kN
     UNIBASES - If UNIBASES is specified, the program expects that a
        unipath basevector is provided, rather than a unipath
        <vecKmerPath>.  See <unibases>.

     SHOW_ALIGNS - when showing longest incorrect unipaths, show their
        alignments to the reference.
     SHOW_PERFECT - If SHOW_PERFECT=n is specified, show the
        alignments of the n largest perfect unipaths.  BADS - show up
        to this many incorrect unipaths

     COPYNUM - if true, show copy number prediction stats
     GC_BINS - if showing the statistics on how well we identify
        copy-number-ploidy unipaths, break down the results by these
        gc-content bin boundaries
     LEN_BINS - if showing the statistics on how well we identify
        copy-number-ploidy unipaths, break down the results by these
        unipath length boundaries

     SHOW_MAP - print out the order of the unipaths along the genome,
        with gaps highlighted unipaths that align imperfectly are
        given a "copy number" of -1, while unipaths that align
        perfectly but wrap around the genome from beginning to end are
        given a "copy number" of -2

     FILTER - filter placements for SHOW_MAP by dominance (i.e. throw
        something out if something longer aligns to the same place).
        Useful if what we're aligning aren't actually unipaths.

     SIZE_UNITS - determines printing of size statistics; options are
        "kmers" or "bases".  For the latter, the total bases printed
        will be exaggerated slightly by overlaps.

     SAVE_FW_PLACEMENTS_TO - if not empty, save fw placement as this
        (relative to run).
     SAVE_RC_PLACEMENTS_TO - if not empty, save rc placement as this
        (relative to run).
     SAVE_ALIGNS - if true, save alignments as a vec of look_aligns.

     ARCHIVE - if given, send output to specified log file, rather
        than cout.
	
   See also <UnipathStats>.

   @file
*/
int main( int argc, char *argv[] )
{
  RunTime( );

  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN); 
  CommandArgument_Int(K);
  CommandArgument_String_OrDefault(GENOME,"genome"); 
  CommandArgument_String_OrDefault(UNIPATHS,"unipaths"); 
  CommandArgument_String_OrDefault(UNIBASES,"");
  CommandArgument_Bool_OrDefault_Doc(UNIPATHS_BIG, False,
				     "whether the unipaths interval database uses big_tagged_rpint "
				     "rather than tagged_rpint" );
  CommandArgument_String_OrDefault(READS,"reads"); 
  CommandArgument_Bool_OrDefault_Doc(SHOW_ALIGNS, False,
				     "show alignments of imperfect unipaths to reference" );
  CommandArgument_Int_OrDefault(SHOW_PERFECT, 0);
  CommandArgument_Int_OrDefault_Doc(BADS, 10,
				    "report on this many incorrect unipaths" );
  CommandArgument_Bool_OrDefault(COPYNUM, False);
  CommandArgument_String_OrDefault(GC_BINS, "{.33,.67}");
  CommandArgument_String_OrDefault(LEN_BINS, "{120}");
  CommandArgument_Bool_OrDefault(SHOW_MAP, False);
  CommandArgument_Bool_OrDefault(FILTER, False);
  CommandArgument_String_OrDefault(SAVE_FW_PLACEMENTS_TO,"");
  CommandArgument_String_OrDefault(SAVE_RC_PLACEMENTS_TO,"");
  CommandArgument_Bool_OrDefault( SAVE_ALIGNS, False );
  CommandArgument_String_OrDefault(SIZE_UNITS, "kmers");
  CommandArgument_Bool_OrDefault(ARCHIVE, False);
  CommandArgument_Bool_OrDefault(VERBOSE, False);
  CommandArgument_Bool_OrDefault(EVAL_GRAPH, False);
  CommandArgument_Int_OrDefault_Doc(PLOIDY, -1,
				    "use PLOIDY for sensible output of binary copy number evaluations");
  CommandArgument_Int_OrDefault_Doc(REF_PLOIDY, 1,
				    "PLOIDY of the supplied reference genome");
  CommandArgument_Int_OrDefault_Doc(MATRIX_LIMIT, 5,
				    "restrict the printout of copy number prediction matrix to the first MATRIX_LIMIT copy numbers");
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  EndCommandArguments;
  
  // Thread control (OMP used by SearchFastb used by GetAlignsFast)
  
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );
  
  
  PRINT( REF_PLOIDY );
  // Define directories.

  String data_dir = PRE + "/" + DATA;
  String run_dir = PRE + "/" + DATA + "/" + RUN;

  // Log stream.

  String log_file = run_dir + "/" + READS + ".";
  log_file += UNIBASES != "" ? UNIBASES : UNIPATHS;
  log_file += ".k" + ToString( K ) + ".eval";

  ofstream arcout;
  if ( ARCHIVE ) arcout.open( log_file.c_str( ) );
  ostream &out = ARCHIVE? arcout : cout;
  if ( ARCHIVE ) {
    cout << "Sending output to " << log_file << "\n" << endl;
    PrintCommandPretty( out );
  }
     
  // Check for ploidy file
    
  if ( PLOIDY == -1 ){
    if ( IsRegularFile( data_dir + "/ploidy") ){
      out << "Reading PLOIDY from the file" << endl;
      PLOIDY = StringOfFile( data_dir + "/ploidy", 1 ).Int( );
      out << "PLOIDY = " << PLOIDY << endl;
    }else{
      InputErr("Set non-negative PLOIDY option or provide \'ploidy\' file in the data directory");
    }
  }
  PRINT( PLOIDY );
  // Load the unipaths.

  out << Date() << " Loading the unipaths/unibases..." << endl;

  vecKmerPath unipaths;
  vecbasevector unibases;
  String strkK = "k" + ToString( K );
  if ( UNIBASES != "" ) {
    String unibaseFile = run_dir + "/" + READS + "." + UNIBASES + "." + strkK;
    if ( !IsRegularFile( unibaseFile ) )
      unibaseFile += ".fastb";
    unibases.ReadAll( unibaseFile );
  }
  if ( (EVAL_GRAPH || UNIBASES == "") && UNIPATHS != "" )
    unipaths.ReadAll( run_dir + "/" + READS + "." + UNIPATHS + ".k" + ToString(K) );
     
  if ( unipaths.empty( ) && unibases.empty( ) ) {
    out << "There are no unipaths or unibases to evaluate.\n";
    exit(0);
  }

  if ( UNIBASES != "" )
  {    vecbasevector genome( data_dir + "/" + GENOME + ".fastb" );
       vec<basevector> missing_genomic_kmers;
       int64_t non_genomic_kmers = 0, N50_unibase = 0;
       if ( K == 64 )
       {    UnibaseSummaryStats<64>( unibases, genome, missing_genomic_kmers,
                 non_genomic_kmers, N50_unibase );    }
       else if ( K == 80 )
       {    UnibaseSummaryStats<80>( unibases, genome, missing_genomic_kmers,
                 non_genomic_kmers, N50_unibase );    }
       else if ( K == 96 )
       {    UnibaseSummaryStats<96>( unibases, genome, missing_genomic_kmers,
                 non_genomic_kmers, N50_unibase );    }
       else if ( K == 128 )
       {    UnibaseSummaryStats<128>( unibases, genome, missing_genomic_kmers,
                 non_genomic_kmers, N50_unibase );    }
       else if ( K == 144 )
       {    UnibaseSummaryStats<144>( unibases, genome, missing_genomic_kmers,
                 non_genomic_kmers, N50_unibase );    }
       else 
       {    cout << "Only implemented for K = 64, 80, 96, 128, 144." << endl;
            cout << "Abort." << endl;
            exit(1);    }
       longlong kmers_in_genome = 0, bases_in_genome = 0;
       for ( size_t t = 0; t < genome.size( ); t++ )
       {    bases_in_genome += genome[t].size( );
            if ( genome[t].isize( ) >= K )
                 kmers_in_genome += genome[t].size( ) - K + 1;    }
       double missing_kmers = 10000.0 * double( missing_genomic_kmers.size( ) )
            / double(kmers_in_genome);
       PerfStat::log( ) << std::fixed << std::setprecision(2) 
            << PerfStat( "missing_kmers", "%% of genomic kmers not in unipaths",
            missing_kmers );
       PerfStat::log( ) << std::fixed << std::setprecision(2) 
            << PerfStat( "unipaths_per_kb", "number of unipaths per kb of genome",
            1000.0 * double( unibases.size( ) ) / double(bases_in_genome) );    }

  // Evaluate Unipath Graph

  if ( EVAL_GRAPH )
    EvalGraph( K, run_dir, READS, UNIPATHS, unipaths, out );
     
#if 0     
  {
    out << Date() << " Loading read paths..." << endl;
    vecKmerPath paths( run_dir + "/" + "reads.paths.k" + ToString(K) );
    vecKmerPath paths_rc( run_dir + "/" + "reads.paths_rc.k" + ToString(K) );
    BREAD2( run_dir + "/reads.pathsdb.k" + ToString(K),
	    vec<tagged_rpint>, pathsdb );
    BREAD2( data_dir + "/genome.unipathsdb_big.k20",
	    vec<big_tagged_rpint>, genome_unipathsdb );
    BREAD2( run_dir + "/reads." + UNIPATHS + "db.k20",
	    vec<tagged_rpint>, unipathsdb );
    vecKmerPath genome_unipaths;
    genome_unipaths.ReadAll( data_dir + "/genome.unipaths.k" + ToString(K) );
       
    out << "Building unipath graph\n";
    digraph A;
    BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb,
				unipaths, unipathsdb, A );
    
    out << Date() << " Loaded genome and read unipaths." << endl;
    EvalUnipaths( unipaths, unipathsdb, genome_unipaths,
		  genome_unipathsdb, A, out );

    out << Date() << " Evaluated genome and read unipaths." << endl;
  }
#endif     
     
  vec<double> gcBins;
  vec<int> lenBins;
  ParseDoubleSet(GC_BINS, gcBins);
  ParseIntSet(LEN_BINS, lenBins);

  // Read in output of UnipathCoverage.

  out << Date() << " Reading in output of UnipathCoverage..." << endl;
  VecPdfEntryVec estimatedCopyNumber;
  String predictedCountFileName =
    run_dir + "/" + READS + "." + UNIPATHS + ".predicted_count.k"
    + ToString(K) ;

  if ( COPYNUM  ) {
    if (!IsRegularFile( predictedCountFileName )) {
      out << "WARNING: disabling copy number prediction stats because the file with " <<
	"copy number predictions, " << predictedCountFileName << ", is not found.\n" << endl;
      COPYNUM = false;
    } else {
      estimatedCopyNumber.ReadAll( predictedCountFileName.c_str() );
      if ( estimatedCopyNumber.size() != (UNIBASES == "" ? unipaths.size() : unibases.size()) ) {
	out << "WARNING: disabling copy number prediction stats - there are " <<
	  unibases.size() << (UNIBASES == "" ? "unipaths" : "unibases")
	    << " in " <<
	  (run_dir + "/" + READS + "." + (UNIBASES=="" ? UNIPATHS : UNIBASES) + ".k" + ToString(K))
	    << " but "
	    << estimatedCopyNumber.size() << " copy number estimates in "
	    << predictedCountFileName << endl;
	COPYNUM = false;
      }
    }
  }
     
  // Local var: copyNumberPloidyStats
  // Statistics on how well we predict copy-number-ploidy-ness of unipaths
  // broken down by unipath length and gc content.
  String propertyYes = "cn==" + ToString(PLOIDY);
  String propertyNo  = "cn!=" + ToString(PLOIDY);
  BinsVec2 < unipath_size_t, gc_t, PredictionStats >
    copyNumberPloidyStats( lenBins, gcBins,
			   PredictionStats( propertyYes, propertyNo ) );

  // Print size statistics.

  vec<int> sizes;
  if ( UNIBASES == "" ) {
    for ( size_t i = 0; i < unipaths.size( ); i++ )
      sizes.push_back( unipaths[i].KmerCount( ) );
  }
  else {
    for ( size_t i = 0; i < unibases.size( ); i++ ) {
      int n = unibases[i].isize( ) - K + 1;
      ForceAssertGe( n, 1 );
      sizes.push_back(n);
    }
  }
  vec<nkmers_t> usizes(sizes);
  Sort(sizes);
  for ( int j = 0; j < sizes.isize( ); j++ )
    sizes[j] += ( SIZE_UNITS == "kmers" ? 0 : K - 1 );
  int zeros = int( ceil( log10( Max(sizes) ) ) );
  vec<int> size_counts( zeros, 0 ), size_counts_total( zeros, 0 );
  for ( int i = 0; i < sizes.isize( ); i++ )
    {    int z = int( ceil( log10( sizes[i] ) ) );
      if ( z > 0 ) --z;
      ++size_counts[z];
      size_counts_total[z] += sizes[i];    }
  out << sizes.size( ) << " unipaths, largest = " << Max(sizes) << " "
      << ( SIZE_UNITS == "kmers" ? "kmers" : "bases" )
      << ", N50 = " << N50( sizes ) << ", total = " << Sum(sizes) << "\n" << endl;
  int bottom = 1, top = 10;
  vec< vec<String> > rows, rows0;
  vec<String> row1, row2, col1;
  row1.push_back( "range of sizes", "number of", "total " + SIZE_UNITS );
  row2.push_back( "", "unipaths" );
  rows.push_back( row1, row2 );
  ostrstream rout;
  for ( int i = 0; i < zeros; i++ )
    {    vec<String> row;
      row.push_back( ToString(bottom), "-", ToString(top) );
      rows0.push_back(row);
      if ( bottom != 1 ) --bottom;
      bottom = 1 + bottom*10, top *= 10;    }
  PrintTabular( rout, rows0, 2, "rrr" );
  rout << ends;
  vec<char> separators;
  separators.push_back( '\n' );
  Tokenize( rout.str( ), separators, col1 );
  for ( int i = 0; i < zeros; i++ )
    {    vec<String> row;
      row.push_back( col1[i], ToString( size_counts[i] ), 
		     ToString( size_counts_total[i] ) );
      rows.push_back(row);    }
  PrintTabular( out, rows, 2, "rrr" );
  out << "\n";

  // Do we know the reference?

  if ( !IsRegularFile( data_dir + "/"+GENOME+".fastb" ) ) {
    out << "No genome to compare to.  Exiting." << endl;
    exit(0);
  }
  if ( !IsRegularFile( data_dir + "/"+GENOME+".lookup" ) ) {
    out << "Lookup table not built for genome.  Exiting." << endl;
    exit(0);
  }

  // Convert unipaths into sequences.

  String KS = ToString(K);
  if ( UNIBASES == "" ) {
    if ( !UNIPATHS_BIG ) {
      KmerBaseBroker kbb( run_dir, K, READS );
      for ( size_t i = 0; i < unipaths.size( ); i++ ) {
	unibases.push_back_reserve( kbb.Seq( unipaths[i] ) );    
      }
    } else {
      KmerBaseBrokerBig kbb( run_dir, K, READS );
      for ( size_t i = 0; i < unipaths.size( ); i++ ) {
	unibases.push_back_reserve( kbb.Seq( unipaths[i] ) );    
      }
    }
  }
  
  // Align the unipath sequences to the reference.

  vecbasevector genome( data_dir + "/" + GENOME + ".fastb" );
  vec<int> mismatches;
  vec<int> indels;
  vec<look_align> aligns1, aligns2, aligns3;
  vec<look_align> to_save;
  vec<look_align> best_aligns;
  vec<genome_placement> fw_placements, rc_placements;
  vec< pair<nkmers_t,int> > len_id;
  vecbitvector cov;
  int nseq = 0;
  
  // create an array of bitvectors parallel to our array of
  // basevectors, to represent what genome bases are covered by
  // unipaths.

  Mimic( genome, cov ); 
  vecbitvector cov100( cov );
  vec<Bool> perfect( unibases.size( ), False );
  
  temp_file uba( "/tmp/UnipathEval_aligns_XXXXXXX" );

  AlignUnibases( K, BADS, PLOIDY, REF_PLOIDY, MATRIX_LIMIT, GENOME, SAVE_ALIGNS,
		 SHOW_MAP, COPYNUM, data_dir, run_dir, gcBins, lenBins, sizes,
		 usizes, estimatedCopyNumber, uba,
		 copyNumberPloidyStats, unibases, aligns1, aligns2, 
		 aligns3, to_save, best_aligns,
		 fw_placements, rc_placements, len_id,
		 cov, cov100, perfect, mismatches,
		 indels, nseq, out );
  
  // Post process aligns3.

  if ( ! len_id.empty( ) ) {
    vec< vec<int> > aligns3_index;
    LoadLookAligns( uba, aligns3, aligns3_index, nseq );
    for ( int seq = 0; seq < nseq; seq++ ) {
      int id = len_id[seq].second;
      int len = len_id[seq].first;
      mismatches[seq] = indels[seq] = len;
      vec<look_align> wrapped_aligns_fw, wrapped_aligns_rc;
      for ( int j = 0; j < aligns3_index[seq].isize( ); j++ ) {
	const look_align& la = aligns3[ aligns3_index[seq][j] ];
	// Ignore improper aligns.
	if ( ( la.a.pos1() != 0 && la.a.pos2() != 0 ) ||
	     ( la.a.Pos1() != (int)la.query_length && 
	       la.a.Pos2() != (int)la.target_length ) ) continue;
	// Defer aligns that wrap unipaths that from end to begin of reference.
	if ( la.a.pos2() == 0 || la.a.Pos2() == (int) la.target_length )
	  ( la.rc1 ? wrapped_aligns_rc : wrapped_aligns_fw ).push_back( la );
	else if ( la.Errors( ) < mismatches[seq] + indels[seq] ) {
	  mismatches[seq] = la.mutations;
	  indels[seq] = la.indels;    
	  if ( SHOW_MAP ) best_aligns[id] = la;
	}    
      }
	 
      for ( int rc = 0; rc < 2; ++rc ) {
	// Find pairs of wrapped aligns that cover the entire unipath.
	vec<look_align>& wrapped_aligns
	  = ( rc==0 ? wrapped_aligns_fw : wrapped_aligns_rc );
	if ( wrapped_aligns.size() < 2 ) continue;
	for ( unsigned int i = 0; i < wrapped_aligns.size()-1; ++i ) {
	  if ( perfect[id] ) break;
	  look_align& aligni = wrapped_aligns[i];
	  for ( unsigned int j = i+1; j < wrapped_aligns.size(); ++j ) {
	    if ( perfect[id] ) break;
	    look_align& alignj = wrapped_aligns[j];
	    // We want aligns from different ends of the same reference.
	    if ( aligni.target_id != alignj.target_id ||
		 aligni.a.pos2() == alignj.a.pos2() ||
		 aligni.a.Pos2() == alignj.a.Pos2() ) 
	      continue;
	    // We want them to meet in the middle of the unipath.
	    if ( aligni.a.Pos1() == alignj.a.pos1() ||
		 aligni.a.pos1() == alignj.a.Pos1() ) {
	      int sumerrs = aligni.Errors() + alignj.Errors();
	      if ( sumerrs < mismatches[i] + indels[i] ) {
		mismatches[seq] = aligni.mutations + alignj.mutations;
		indels[seq] = aligni.indels + alignj.indels;
	      }
	      if ( sumerrs == 0 ) {
		perfect[id] = true;
		for ( int a = 0; a < 2; ++a ) {
		  look_align& la = ( a == 0 ? aligni : alignj );
		  if ( SAVE_ALIGNS ) to_save.push_back( la );
		  for ( int k = la.pos2(); k < la.Pos2( ); k++ )
		    cov[la.target_id].Set( k, True );
		  if ( usizes[id] > 100 )
		    for ( int k = la.pos2(); k < la.Pos2( ); k++ )
		      cov100[la.target_id].Set( k, True );
		  if ( SHOW_MAP ) {
		    genome_placement pl( id, la.extent1()-K+1, 
					 la.target_id, la.pos2(), 
					 la.Pos2(), la.rc1, -2 );
		    if ( la.rc1 )
		      rc_placements.push_back( pl );
		    else
		      fw_placements.push_back( pl );
		  }
		}
	      } // sumerrs == 0
	    } // meet in middle
	  } // j loop
	} // i loop
      } // rc loop
    } // seq loop
  } // len_id != empty

  // Report.

  out << "base coverage by perfect unipaths = "
      << setprecision(4) << 100.0 * Coverage(cov) << "%\n";
  out << "base coverage by perfect unipaths of size > 100 = "
      << setprecision(4) << 100.0 * Coverage(cov100) << "%\n";
  int total100 = 0, good100 = 0;
  size_t good = 0;
  for ( size_t i = 0; i < unibases.size( ); i++ ) {
    if ( usizes[i] > 100 ) {
      ++total100;
      if ( perfect[i] ) ++good100;
    }
    if ( perfect[i] ) ++good;
  }
  out << PERCENT_RATIO( 4, good, unibases.size( ) )
      << " of all unipaths are perfect\n";
  out << PERCENT_RATIO( 4, good100, total100 )
      << " of unipaths of size > 100 are perfect\n\n";
  if ( nseq == 0 ) out << "All unipaths match reference perfectly.\n";
  else out << "longest incorrect unipaths:\n";

  for ( int i = 0; i < nseq; i++ ) {
    int id = len_id[i].second, len = len_id[i].first;
    out << "[" << i << "] unipath " << id << " (length=" << len << ")";
    if ( perfect[id] ) {
      out << " - actually perfect, but connects the ends of the reference\n";
      continue;
    }
    out << " - best placement on reference has ";
    if ( indels[i] == 0 ) out << mismatches[i] << " mismatches\n";
    else if ( mismatches[i] == 0 ) out << indels[i] << " indels\n";
    else {
      out << mismatches[i] << " mismatches and " << indels[i] 
	  << " indels\n";
    }
  }
  out << "\n";
  if (SHOW_ALIGNS) {
    String line;
    Ifstream( ain, uba );
    while(1) {
      getline( ain, line );
      if ( !ain ) break;
      if ( line.Contains( "QUERY", 0 ) ) continue;
      out << line << "\n";
    }
  }
  if ( SHOW_PERFECT > 0 ) {
    out << "\nalignments of longest perfect unipaths:\n\n";
    vec< pair<nkmers_t,align_id_t> > len_ind;
    for ( align_id_t i = 0; i < aligns1.isize( ); i++ )
      len_ind.push( usizes[ aligns1[i].query_id ], i );
    ReverseSort(len_ind);
    for ( int i = 0; i < Min( SHOW_PERFECT, len_ind.isize( ) ); i++ ) {
      const look_align& la = aligns1[ len_ind[i].second ];
      out << "unipath " << la.query_id 
	  << " (length=" << usizes[la.query_id]
	  << ") aligns to " << la.target_id << "." << la.pos2( ) << "-"
	  << la.Pos2( ) << "\n\n";
    }
  }

  if ( COPYNUM && estimatedCopyNumber.size() > 0 ) {
    nkmers_t maxUnipathSize = Max( usizes );
       
    copyNumberPloidyStats.SetOuterBounds(0, maxUnipathSize, 0, 1.0);
    for (bin_id_t lenBin = 0; lenBin <= lenBins.size(); lenBin++) {
      for (bin_id_t gcBin = 0; gcBin <= gcBins.size(); gcBin++) {
	unipath_size_t lenLo, lenHi;
	gc_t gcLo, gcHi;
	copyNumberPloidyStats.GetBinBounds( lenBin, gcBin, lenLo,
					    lenHi, gcLo, gcHi );
	if ( VERBOSE || ARCHIVE ){
	  out << endl << " copy-number-ploidy unipath calls for length in ["
	      << lenLo << "," << lenHi << "), gc in ["
	      << gcLo << "," << gcHi << "):\n\n";
	  
	  copyNumberPloidyStats[lenBin][gcBin].Print( out );
	}
      }
    }
  }

  if ( SHOW_MAP || SAVE_ALIGNS ) {
    if ( SHOW_MAP ) out << "Map of unipaths to genome:" << endl;
       
    vec<copy_num_t> copyNumber( unibases.size(), 0 );
    for ( int rc = 0; rc < 2; ++rc ) {
      vec<genome_placement>& placements
	= ( rc ? rc_placements : fw_placements );
      for ( int i = 0; i < placements.isize(); ++i )
	copyNumber[ placements[i].GetReadId() ]++;
    }
    for ( int rc = 0; rc < 2; ++rc ) {
      vec<genome_placement>& placements
	= ( rc ? rc_placements : fw_placements );
      for ( int i = 0; i < placements.isize(); ++i )
	if ( placements[i].GetCopyNumber() == 0 )
	  placements[i].SetCopyNumber( copyNumber[placements[i].GetReadId()] );
    }
       
    for ( size_t i = 0; i < unibases.size(); ++i ) {
      look_align& la = best_aligns[i];
      if ( la.target_length > 0 ) {
	genome_placement pl( i, unibases[i].size()-K+1, la.target_id,
			     la.pos2(), la.Pos2(), la.rc1, -1 );
	if ( SAVE_ALIGNS ) {
	  // The id of unipath to be saved is i.
	  look_align saved_la = la;
	  saved_la.query_id = i; 
	  to_save.push_back( saved_la ); 
	}
	if ( la.rc1 ) rc_placements.push_back( pl );
	else fw_placements.push_back( pl );
      }
    }

    for ( int rc = 0; rc < 2; ++rc ) {
      vec<genome_placement>& placements
	= ( rc ? rc_placements : fw_placements );
	 
      sort( placements.begin(), placements.end() );
	 
      if ( FILTER ) {
	vec<Bool> filtered( placements.size(), False );
	list<genome_placement> current;
	for ( int i = 0; i < placements.isize(); ++i ) {
	  list<genome_placement>::iterator c = current.begin();
	  while ( c != current.end() )
	    if ( c->GetEndOnGenome() <= placements[i].GetStartOnGenome() ) {
	      list<genome_placement>::iterator dead = c++;
	      current.erase(dead);
	    }
	    else {
	      if ( c->GetStartOnGenome() < placements[i].GetStartOnGenome() &&
		   c->GetEndOnGenome() >= placements[i].GetEndOnGenome() ||
		   c->GetStartOnGenome() <= placements[i].GetStartOnGenome() &&
		   c->GetEndOnGenome() > placements[i].GetEndOnGenome() ) {
		filtered[i] = True;
		break;
	      }
	      c++;
	    }
	  if ( ! filtered[i] )
	    current.push_back( placements[i] );
	}
           
	EraseIf( placements, filtered );
      }

      if ( SHOW_MAP ) {
	out << ( rc ? "RC MAP" : "FW MAP" ) << endl;
	   
	if ( placements.empty() ) {
	  out << "No alignments found!" << endl;
	} else {
	  out << placements.front();
	  for ( int i = 1; i < placements.isize(); ++i ) {
	    nbases_t gapInBases
	      = placements[i].GetStartOnGenome()
	      - placements[i-1].GetEndOnGenome();
	    nkmers_t gapInKmers = gapInBases + K-1;
	    if ( gapInKmers > 0 )
	      out << "  GAP of " << gapInKmers << " kmers" << endl;
	    out << placements[i];
	  }
	}
      }
    }
  }
     
  if ( ! SAVE_FW_PLACEMENTS_TO.empty() )
    BinaryWrite3( run_dir + "/" + SAVE_FW_PLACEMENTS_TO, fw_placements );
     
  if ( ! SAVE_RC_PLACEMENTS_TO.empty() )
    BinaryWrite3( run_dir + "/" + SAVE_RC_PLACEMENTS_TO, rc_placements );

  if ( SAVE_ALIGNS ) {
    out << "\n" << Date( ) << ": sorting and saving alignments\n" << endl;
       
    sort( to_save.begin( ), to_save.end( ) );
       
    String alout_file = log_file.Before( ".eval" ) + ".qlt";
    ofstream alout( alout_file.c_str( ) );
    for (int ii=0; ii<to_save.isize( ); ii++) {
      const look_align &al = to_save[ii];
      const basevector *qbases = &( unibases[al.QueryId( )] );
      const basevector *tbases = &( genome[al.TargetId( ) ] );
      al.PrintParseable( alout, qbases, tbases );
    }
    alout.close( );
  }
  
  out << Date() << ": Done with UnipathEval!" << endl;
}



/***************************************************************************
 *
 * Functions start here.
 *
 ***************************************************************************/



/**
 * Function: EvalUnipaths
 *
 * Evaluate unipaths by taking each genomic unipath, and measuring how
 * much it got broken up due to missing kmer adjacencies.  As a
 * variation, we may measure this for just the "important" unipaths,
 * for example long-enough CN1 unipaths.
 */
void EvalUnipaths(  const vecKmerPath& unipaths,
		    const vec<tagged_rpint>& unipathsdb,
		    const vecKmerPath& genome_unipaths,
		    const vec<big_tagged_rpint>& genome_unipathsdb,
		    const digraph& unipathAdjGraph,
		    ostream &out )
{
  out << " there are " << unipaths.size() << " read unipaths, with " << unipathsdb.size() << " intervals." << endl;
  out << " there are " << genome_unipaths.size() << " genome unipaths, with " << genome_unipathsdb.size() << " intervals." << endl;

  CheckUnipathSanity( unipaths, unipathsdb );
  CheckUnipathSanity( genome_unipaths, genome_unipathsdb );

  VecLongVec genome_unipath_interval_coverers( genome_unipathsdb.size() );
  path_interval_id_t unipathsdb_idx = 0;
  int maxCoverers = 0;
  int minCoverers = 0;
  for ( path_interval_id_t genome_unipathsdb_idx = 0;
	genome_unipathsdb_idx < genome_unipathsdb.isize() &&
	  unipathsdb_idx < unipathsdb.isize(); genome_unipathsdb_idx++ ) {
    while ( unipathsdb_idx < unipathsdb.isize() &&
	    unipathsdb[ unipathsdb_idx ].Stop() < genome_unipathsdb[ genome_unipathsdb_idx ].Start() )
      unipathsdb_idx++;
    while ( unipathsdb_idx < unipathsdb.isize() &&
	    IntervalsOverlap( unipathsdb[ unipathsdb_idx ].Start(), unipathsdb[ unipathsdb_idx ].Stop(),
			      genome_unipathsdb[ genome_unipathsdb_idx ].Start(), genome_unipathsdb[ genome_unipathsdb_idx ].Stop() ) ) {
      ForceAssert( genome_unipath_interval_coverers[ genome_unipathsdb_idx ].empty()  ||
		   genome_unipath_interval_coverers[ genome_unipathsdb_idx ].back() < unipathsdb_idx );
      genome_unipath_interval_coverers[ genome_unipathsdb_idx ].push_back( unipathsdb_idx++ );
    }

  }

  // For each genome unipath, make a list of indices of its intervals.
  VecLongVec genome_unipath_interval_ids( genome_unipaths.size() );
  for ( path_interval_id_t genome_unipathsdb_idx = 0; genome_unipathsdb_idx < genome_unipathsdb.isize(); genome_unipathsdb_idx++ ) {
    unipath_id_t genome_unipath_id = genome_unipathsdb[ genome_unipathsdb_idx ].ReadId();
    LongVec& this_genome_unipath_interval_ids = genome_unipath_interval_ids[ genome_unipath_id ];
    if ( this_genome_unipath_interval_ids.empty() )
      this_genome_unipath_interval_ids.resize( genome_unipaths[ genome_unipath_id ].NSegments() );
    this_genome_unipath_interval_ids[ genome_unipathsdb[ genome_unipathsdb_idx ].PathPos() ] = genome_unipathsdb_idx;
  }

  nkmers_t diffs = 0;
  
  // Now go through each genome unipath
  for ( size_t genome_unipath_id = 0; genome_unipath_id < genome_unipaths.size(); genome_unipath_id++ ) {

    unipath_id_t cur_unipath = -1;
    nkmers_t cur_unipath_overlap_len = 0;

    ForceAssert( genome_unipath_interval_ids[ genome_unipath_id ].size() ==
                static_cast<unsigned>(genome_unipaths[ genome_unipath_id ].NSegments()) );

    vec< nkmers_t > overlap_lens;
    
    for ( unsigned int gseg = 0; gseg < genome_unipath_interval_ids[ genome_unipath_id ].size(); gseg++ ) {
      path_interval_id_t genome_unipathsdb_idx = genome_unipath_interval_ids[ genome_unipath_id][ gseg ];
      for ( unsigned int j = 0; j < genome_unipath_interval_coverers[ genome_unipathsdb_idx ].size(); j++ ) {
	path_interval_id_t unipathsdb_idx = genome_unipath_interval_coverers[ genome_unipathsdb_idx ][j];

	ForceAssert( IntervalsOverlap( unipathsdb[ unipathsdb_idx ].Start(), unipathsdb[ unipathsdb_idx].Stop(),
				       genome_unipaths[ genome_unipath_id ].Start( gseg ), genome_unipaths[ genome_unipath_id ].Stop( gseg ) ) );

	nkmers_t overlapLen =
	  Min( unipathsdb[ unipathsdb_idx].Stop(), genome_unipaths[ genome_unipath_id ].Stop( gseg ) ) -
	  Max( unipathsdb[ unipathsdb_idx].Start(), genome_unipaths[ genome_unipath_id ].Start( gseg ) )
	  + 1;
	ForceAssert( overlapLen > 0 );

	unipath_id_t this_unipath = unipathsdb[ unipathsdb_idx ].ReadId();
	if ( this_unipath == cur_unipath )
	  cur_unipath_overlap_len += overlapLen;
	else if ( cur_unipath >= 0  && unipathAdjGraph.HasEdge( cur_unipath, this_unipath ) ) {
	  cur_unipath_overlap_len += overlapLen;
	  cur_unipath = this_unipath;
	} else {
	  if ( cur_unipath_overlap_len > 0 )
	    overlap_lens.push_back( cur_unipath_overlap_len );
	  cur_unipath = unipathsdb[ unipathsdb_idx ].ReadId();
	  cur_unipath_overlap_len = overlapLen;
	}
      }  // for each coverer of this genome unipath interval
    }  // for each segment of this genome unipath
    
    if ( cur_unipath_overlap_len > 0 )
      overlap_lens.push_back( cur_unipath_overlap_len );

    Sort( overlap_lens );
    nkmers_t diffHere = genome_unipaths[ genome_unipath_id ].KmerCount() -
      ( overlap_lens.empty() ? 0 : N50( overlap_lens ) );
    ForceAssert( diffHere >= 0 );
    diffs += diffHere;
  }  // for each genome unipath

  PRINT_TO(out, diffs);
}



/**
 * Function: EvalGraph
 */
void EvalGraph ( const int K,
		 const String &run_dir,
		 const String &READS,
		 const String &UNIPATHS,
		 const vecKmerPath &unipaths,
		 ostream &out )
{

  out << Date() << " Loading read paths and database files." << endl;
  BREAD2( run_dir + "/" + READS + "." + UNIPATHS + "db.k" + ToString(K), vec<tagged_rpint>, unipathsdb );
  vecKmerPath paths( run_dir + "/" + READS + ".paths.k" + ToString(K) );
  vecKmerPath paths_rc( run_dir + "/" + READS + ".paths_rc.k" + ToString(K) );    
  BREAD2( run_dir + "/" + READS + ".pathsdb.k" + ToString(K), vec<tagged_rpint>, pathsdb );
            
  out << Date() << " Building unipath graph" << endl;
  digraph A;
  BuildUnipathAdjacencyGraph( paths, paths_rc, pathsdb, unipaths, unipathsdb, A );

  // Determine unipath graph components
  vec< vec<int> > comps;
  A.Components(comps);
    
  pair<int, int> disconnected, connected;
  vec<longlong> compKmerCount( comps.isize(), 0 ); // component size in kmers
  vec<int> compUnipathCount(comps.isize());   // number of unipaths in component
  for (int i = 0; i < comps.isize(); i++) {
    compUnipathCount[i] = comps[i].isize();
    int kmerCount = 0;
    for (int j = 0; j < compUnipathCount[i]; j++){
      kmerCount += unipaths[comps[i][j]].KmerCount();
    }
    compKmerCount[i] = kmerCount;
    if (compUnipathCount[i] == 1) {
      disconnected.first++;
      disconnected.second += kmerCount;
    } else {
      connected.first += compUnipathCount[i];
      connected.second += kmerCount;
    }
  }

  longlong totalKmerCount = disconnected.second + connected.second;
  out << "Unipath graph consists of " << comps.isize() << " components (" << totalKmerCount << " kmers)" << endl;
  out << "Disconnected unipaths : " << disconnected.first << " (" << disconnected.second << " kmers)" << endl;
  out << "Connected unipaths   : " << connected.first << " (" << connected.second << " kmers)" << endl;
       
  ReverseSortSync( compKmerCount, compUnipathCount );
  out << "20 largest unipath graph component (in kmers):" << endl;
  for ( int i = 0; i < Min(compKmerCount.isize(), 20); i++ ) {
    double percSize = 100.0 * (double) compKmerCount[i] / (double) totalKmerCount;
    out << "\t" << compUnipathCount[i] << " unipaths,  " << compKmerCount[i] << " kmers,  " << percSize << " % of all kmers" << endl;
  }

}



/**
 * Function: EvalCopyNumber
 */
void EvalCopyNumber( const int PLOIDY,
		     const int MATRIX_LIMIT,
		     const vecbvec &unibases,
		     const VecPdfEntryVec &estimatedCopyNumber,
		     const vec<copy_num_t> &trueCopyNumber,
		     const vec<double> &gcBins,
		     const vec<int> &lenBins,
		     const vec<int> &sizes,
		     BinsVec2 < unipath_size_t, gc_t, PredictionStats > &copyNumberPloidyStats,
		     ostream &out )
{
  out << "\n            ------ COPY NUMBER EVALUATION ------" << endl;
  out << " NOTE! copy number evaluation is based on perfect lookup for now.\n" << endl;
  vec< pair< nbases_t, prob_t > > actualYesPredictedNo, actualNoPredictedYes;
  vec<copy_num_t> predCopyNumber( unibases.size(), 0 );
  int maxTrueCopyNumber = 0;
  int maxPredCopyNumber = 0;
  for ( VecPdfEntryVec::size_type unipathId = 0;
          unipathId < estimatedCopyNumber.size();
          unipathId++ ) {
    copy_num_t mostProbableCopyNumber = -1;
    prob_t mostProbableCopyNumberProb = -1.0;
    for (PdfEntryVec::size_type copyNumberProbEntry = 0;
            copyNumberProbEntry < estimatedCopyNumber[unipathId].size();
            copyNumberProbEntry++) {
      const pdf_entry& copyNumberProb = estimatedCopyNumber[unipathId][copyNumberProbEntry];
      copy_num_t thisCopyNumber = copyNumberProb.NumCopies();
      prob_t thisCopyNumberProb = copyNumberProb.Prob();
      if ( thisCopyNumberProb > mostProbableCopyNumberProb ) {
	mostProbableCopyNumber = thisCopyNumber;
	mostProbableCopyNumberProb = thisCopyNumberProb;
      }
    }
    
    double gcPercent = GcPercent( unibases[unipathId] ) / 100.;
    copyNumberPloidyStats(unibases[unipathId].size(), gcPercent)(mostProbableCopyNumber == PLOIDY, trueCopyNumber[unipathId] == PLOIDY)++;
	 
    if ( mostProbableCopyNumber == PLOIDY && 
	 trueCopyNumber[unipathId] != PLOIDY )
      actualNoPredictedYes.push( unibases[ unipathId ].size(), mostProbableCopyNumberProb );
	 
    if ( mostProbableCopyNumber != PLOIDY && 
	 trueCopyNumber[unipathId] == PLOIDY )
      actualYesPredictedNo.push( unibases[ unipathId ].size(), mostProbableCopyNumberProb );
	 
    if ( mostProbableCopyNumber == PLOIDY && 
	 trueCopyNumber[unipathId] > PLOIDY ){
      out << "WARNING! ";
      PRINT5_TO( out,unipathId, unibases[unipathId].size(), 
		 trueCopyNumber[unipathId], mostProbableCopyNumber,
		 mostProbableCopyNumberProb );
    }

    predCopyNumber[ unipathId ] = mostProbableCopyNumber;
    if ( mostProbableCopyNumber > maxPredCopyNumber ) 
      maxPredCopyNumber = mostProbableCopyNumber;
    if ( trueCopyNumber[ unipathId ] > maxTrueCopyNumber )
      maxTrueCopyNumber = trueCopyNumber[ unipathId ];

  }

  ReverseSort( actualNoPredictedYes );
  ReverseSort( actualYesPredictedNo );       

  int numToShow = min(actualYesPredictedNo.isize(), 10);
  if ( numToShow == 0 )
    out << "\nNo unipaths with actual CN-PLOIDY copy numbers were mispredicted. " << endl;
  else{
    out << "\nActual CN-PLOIDY, predicted not CN-PLOIDY: showing first " << numToShow 
	<< " out of " << actualYesPredictedNo.isize() << endl;
    for ( int i = 0; i < numToShow; i++ )
      out << " len=" << actualYesPredictedNo[ i ].first << " conf=" << actualYesPredictedNo[ i ].second << endl;
  }
       
  numToShow = min(actualNoPredictedYes.isize(), 10);
  if ( numToShow == 0 )
    out << "\nNo unipaths with actual non CN-PLOIDY copy numbers were mispredicted as having CN-PLOIDY. " << endl;
  else{
    out << "\nActual not CN-PLOIDY, predicted CN-PLOIDY: showing first " << numToShow
	<< " out of " << actualNoPredictedYes.isize() << endl;
    for ( int i = 0; i < numToShow; i++ ) {
      PRINT_TO( out, actualNoPredictedYes[i].first );
      out << " len=" << actualNoPredictedYes[ i ].first << " conf=" << actualNoPredictedYes[ i ].second << endl;
    }
  }
  
  Bool PrintMatrix = True;

  // add minimum and maximum bounds to bins
  vec<double> gcBins2 = gcBins;
  vec<int> lenBins2 = lenBins;
  gcBins2.push_back( 0.0, 1.001 );
  UniqueSort( gcBins2 );			
  gcBins2.Println(out);
  lenBins2.push_back( 0, Max(sizes) + 1 );
  UniqueSort( lenBins2 );			
  lenBins2.Println(out);
       
  PRINT2( maxTrueCopyNumber, maxPredCopyNumber );
  int maxCopyNumber = max( maxTrueCopyNumber, maxPredCopyNumber );
  
  
  // print information about number of predictions
  vec<int> predCount( maxCopyNumber +1, 0 );
  for ( int u = 0; u < predCopyNumber.isize(); u++ )
    predCount.at( predCopyNumber[u] )++;
  vec<int> trueCount( maxCopyNumber +1, 0 );
  for ( int u = 0; u < trueCopyNumber.isize(); u++ )
    trueCount.at( trueCopyNumber[u] )++;
  out << "\n ---------- numbers of true and predicted copy numbers (cn) -------\n";
  out << "cn\t#true\t#pred" << endl;
  for ( int cn = 0; cn <= maxCopyNumber; cn++ )
    if ( trueCount.at( cn ) > 0 || predCount.at( cn ) > 0 )
      out << cn << "\t" << trueCount[ cn ] << "\t" << predCount[ cn ] << endl;


  out << "\n         -----------------------\n\n";
  out << "prediction matrices (column number indicates alignment based cn, row number indicates predicted cn):" << endl;  
  out << "note: we show only first " << MATRIX_LIMIT << " copy numbers in a matrix\n";
  for ( int iGcBin = 1; iGcBin < gcBins2.isize(); iGcBin++ ){
    double gcInf = gcBins2[ iGcBin -1 ], gcSup = gcBins2[ iGcBin ];
    for ( int iLenBin = 1; iLenBin < lenBins2.isize(); iLenBin++ ){
      double lenInf = lenBins2[ iLenBin -1 ], lenSup = lenBins2[ iLenBin ];
	   
      vec< vec< long > > predMat( maxCopyNumber +1, vec<long> ( maxCopyNumber +1, 0 ) );
      for ( VecPdfEntryVec::size_type unipathId = 0;
              unipathId < estimatedCopyNumber.size();
              unipathId++ ) {
	     
	double uGcContent = GcPercent( unibases[unipathId] ) / 100.0;
	if ( uGcContent < gcInf || uGcContent >= gcSup )  
	  continue;
	int ulen          = unibases[unipathId].isize();
	if ( ulen < lenInf || ulen >= lenSup )  
	  continue;
	     
	     
	int trueCN = trueCopyNumber.at( unipathId );
	int predCN = predCopyNumber.at( unipathId );
	predMat.at( predCN ).at( trueCN ) += 1;
      }
	   
      out << "\n prediction matrix for gc-content [" << gcInf << "," << gcSup << ")  lengths [" << lenInf << "," << lenSup << "):\n  \t";
      for ( int i = 0; i < min( maxCopyNumber +1, MATRIX_LIMIT ); i++ )
	out << i << "\t";
      out << "\n";
      for ( int predCN = 0; predCN < min(predMat.isize(), MATRIX_LIMIT); predCN++ ){
	out << predCN << " \t";
	for ( int trueCN = 0; trueCN < min( predMat[predCN].isize(), MATRIX_LIMIT ); trueCN++ ){
	  out << predMat.at(predCN).at(trueCN) << "\t";
	}
	out << "\n";
      }
	   
    }
  }

}



/**
 * Function: AlignUnibases
 */
void AlignUnibases( const int &K,
		    const int &BADS,
		    const int &PLOIDY,
		    const int &REF_PLOIDY,
		    const int &MATRIX_LIMIT,
		    const String &GENOME,
		    const bool &SAVE_ALIGNS,
		    const bool &SHOW_MAP,
		    const bool &COPYNUM,
		    const String &data_dir,
		    const String &run_dir,
		    const vec<double> &gcBins,
		    const vec<int> &lenBins,
		    const vec<int> &sizes,
		    const vec<nkmers_t> &usizes,
		    const VecPdfEntryVec &estimatedCopyNumber,
		    temp_file &uba,
		    bin2PredStat &copyNumberPloidyStats,
		    vecbvec &unibases,
		    vec<look_align> &aligns1,
		    vec<look_align> &aligns2,
		    vec<look_align> &aligns3,
		    vec<look_align> &to_save,
		    vec<look_align> &best_aligns,
		    vec<genome_placement> &fw_placements,
		    vec<genome_placement> &rc_placements,
		    vec< pair<nkmers_t,int> > &len_id,
		    vecbitvector &cov,
		    vecbitvector &cov100,
		    vec<Bool> &perfect,
		    vec<int> &mismatches,
		    vec<int> &indels,
		    int &nseq,
		    ostream &out )
{
  // Temp dir
  Mkdir777(run_dir + "/tmp");

  vec<copy_num_t> trueCopyNumber( unibases.size( ), 0 );

  PerfectLookup( 12, unibases, data_dir + "/" + GENOME + ".lookup",
		 aligns1, FW_OR_RC );

  for ( int i = 0; i < aligns1.isize( ); i++ ) {
    const look_align& la = aligns1[i];
    unipath_id_t unipathId = la.query_id;
    genome_part_id_t genomePartId = la.target_id;
    if ( SAVE_ALIGNS ) to_save.push_back( la );
    if ( la.FullLength() ) {
      perfect[unipathId] = True;
      if ( SHOW_MAP ) {
	genome_placement pl( unipathId, la.extent1()-K+1, genomePartId,
			     la.pos2(), la.Pos2(), la.rc1, 0 );
	if ( la.rc1 )
	  rc_placements.push_back( pl );
	else
	  fw_placements.push_back( pl );
      }
      if ( COPYNUM ){
	if ( REF_PLOIDY == PLOIDY )
	  trueCopyNumber[unipathId]++;
	else if ( REF_PLOIDY == 1 )
	  trueCopyNumber[unipathId] += PLOIDY;
	else
	  FatalErr("REF_PLOIDY=" + ToString(REF_PLOIDY) + " has unexpected value");
	
      }
    }
    for ( genome_part_pos_t j = la.pos2( ); j < la.Pos2( ); j++ ) {
      cov[genomePartId].Set( j, True );
      if ( usizes[unipathId] > 100 )
	cov100[genomePartId].Set( j, True );
    }
  }
     
  // Evaluate Copy Number

  if ( COPYNUM )
    EvalCopyNumber( PLOIDY, MATRIX_LIMIT, unibases, estimatedCopyNumber,
		    trueCopyNumber, gcBins, lenBins, sizes,
		    copyNumberPloidyStats, out );
     
  // kill perfectly aligned unipaths so that their alignment to the
  // reference will not get caught in subsequent alignment steps.
  // this way we prevent the same alignment from being counted
  // multiple times in various statistics.

  for ( size_t ii=0; ii<unibases.size(); ii++ )
    if ( perfect[ii] )
      unibases[ii].resize(0);
  unibases.WriteAll( run_dir + "/unibases.UnipathEval.fastb" );

  out << Date() << ": Calling GetAlignsFast to align "
      << unibases.size() << " unipaths to the reference..." << endl;
  
  best_aligns.resize( unibases.size() );
  
  // Align the unibases to reference using GetAlignsFast.
  GetAlignsFast( K, run_dir + "/unibases.UnipathEval.fastb",
		 data_dir + "/" + GENOME + ".lookup",
		 run_dir + "/unibases.UnipathEval.aligns",
		 aligns2, false , run_dir + "/tmp");
  LoadLookAligns( run_dir + "/unibases.UnipathEval.aligns", aligns2 );

  out << Date() << ": GetAlignsFast done." << endl;
  for ( size_t i = 0; i < unibases.size( ); i++ )
    if ( !perfect[i] )
      len_id.push_back( make_pair( usizes[i], i ) ); 
  
  ReverseSort(len_id);
  nseq = Min( BADS, len_id.isize( ) );
  mismatches.resize(nseq, 0);
  indels.resize(nseq, 0);
  if ( ! len_id.empty( ) ) {
    temp_file ubf( run_dir + "/tmp/UnipathEval_fasta_XXXXXXX" );
    {
      Ofstream( ubfout, ubf );
      for ( int i = 0; i < nseq; i++ )
	unibases[ len_id[i].second ].Print( ubfout, i );
    }

// MakeDepend: dependency QueryLookupTable
    SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 WE=100 SEQS=" + ubf
	    + " L=" + data_dir + "/" + GENOME + ".lookup"
	    " TMP_DIR=" + run_dir + "/tmp"
	    " PARSEABLE=True VISUAL=True NH=True QUIET=True > " 
	    + uba );
  } // len_id != empty
  
}
