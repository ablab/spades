///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/**
   Program: UnibaseCopyNumber3

   Determines copy number of each unibase using k-mer counts (experimental).
   Adjusts for GC bias.

   INPUT FILES (in run dir):
     <READS_EC>.fastb
     <READS>.<UNIBASES>.k<K>

   OUTPUT FILES:
     reads.<UNIPATHS>.predicted_count.k<K>
   
   CACHED FILES:
     <READS_EC>.<K>merParcels
     <READS>.<UNIBASES>.<K>merParcels
   
   @file
*/

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "paths/PdfEntry.h"
#include "paths/UnibaseUtils.h"
#include "paths/UnibaseCopyNumber3Core.h"

#include "paths/UnibaseCopyNumber3LowerCN.h"



// Check if the gap size from u1 to u2 is the same as from u2* to u1*.
void CheckGaps(
     const vec<int> & to_rc,
     const vec<ulink_with_uids>& condensed_links
    )
{
  map< pair<int,int>, int > gaps;
  for( int i =0; i < condensed_links.isize(); i++ ){
    const ulink_with_uids& link = condensed_links[i];
    int u1 = link.u1;
    int u2 = link.u2;
    gaps[ make_pair(u1,u2) ] = link.sep;
  }
  for ( map< pair<int,int>, int >::iterator it = gaps.begin();
      it != gaps.end(); it++ ) {
    if ( gaps[ make_pair( to_rc[it->first.second], to_rc[it->first.first] ) ] != it->second ) 
      cout << "warning not matching gap " << it->first.first << "_" << it->first.second << 
	"  to  " << to_rc[it->first.second] << "_" <<  to_rc[it->first.first]  << endl;
  }
}

int main( int argc, char** argv ) {
  
  RunTime();
  
  BeginCommandArguments;
  CommandArgument_String(PRE);
  CommandArgument_String(DATA);
  CommandArgument_String(RUN);
  CommandArgument_Int(K);
  CommandArgument_Int_OrDefault( PLOIDY, -1 );
  
  // Infixes for input/output file names.
  CommandArgument_String_OrDefault( READS, "reads" );
  CommandArgument_StringSet_OrDefault( READS_EC, "{reads_ec}" );
  CommandArgument_String_OrDefault( UNIPATHS, "unipaths" );
  CommandArgument_String_OrDefault( UNIBASES, "unibases" );
  
  // Heuristics for bias computation.
  CommandArgument_Bool_OrDefault(CORRECT_BIAS, CORRECT_BIAS_DEFAULT);
  CommandArgument_Double_OrDefault(THRESH, THRESH_DEFAULT);
  CommandArgument_Double_OrDefault(ERR_RATE, ERR_RATE_DEFAULT);
  CommandArgument_Bool_OrDefault(LOWER, False);
  CommandArgument_Int_OrDefault(LOWER_VERBOSITY, 0);

  // Experimental features
  CommandArgument_Bool_OrDefault(EXP_GAP_REMODEL, False);
  CommandArgument_Bool_OrDefault(GAP_STATS_REMODEL, False);
  CommandArgument_Bool_OrDefault(WRITE_GAP, False);
  CommandArgument_String_OrDefault(WRITE_GAP_SUFFIX, "");
  
  // Runtime control.

  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault(WRITE, True);
  CommandArgument_Bool_OrDefault(VERBOSE, False);

  CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
  CommandArgument_Int_OrDefault(MAX_PLACEMENTS, 50);
  
  EndCommandArguments;

  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);
  omp_set_num_threads( NUM_THREADS );

  // Define directories.
  
  String data_dir = PRE + "/" + DATA;
  String run_dir = data_dir + "/" + RUN;

  String merged_head = run_dir + "/" + READS + ".merged_reads";
  
  // Make sure ploidy file (or values) are available
  if ( PLOIDY <= 0 ){
    if ( IsRegularFile( run_dir + "/ploidy" ) )
      PLOIDY = StringOfFile( run_dir + "/ploidy", 1 ).Int( );
    else 
      FatalErr("Require genome ploidy. Please supply either a ploidy "
	       "file in the data, or a value for the PLOIDY option");
    ForceAssertGt( PLOIDY, 0 );
  }

  // Input files
  String unibases_head = run_dir + "/" + READS + "." + UNIBASES;
  String uni_file      = unibases_head + ".k" + ToString(K);
  if ( ! IsRegularFile( uni_file ) )
    InputErr( uni_file + " not found");


  for ( int i = 0; i < READS_EC.isize(); i++) 
    if ( ! IsRegularFile( run_dir + "/" +  READS_EC[i] + ".fastb") )
	 InputErr( run_dir + "/" + READS_EC[i] + ".fastb not found");
  
  String jreads_file = run_dir + "/" + JUMP_READS + ".fastb";
  if ( LOWER && ! IsRegularFile( jreads_file ) )
       InputErr( jreads_file + " not found");

  String nonjumps_file = run_dir + "/" + JUMP_READS + ".outies";
  if ( LOWER && ! IsRegularFile( nonjumps_file ) )
       InputErr( nonjumps_file + " not found");
       
  

  // Load unibases.
  cout << Date() << ": Loading unibases" << endl;
  vecbasevector unibases(uni_file);
  
  
  // Load basevectors
  vecbasevector reads;
  for (int i = 0; i < READS_EC.isize(); i++) {
    cout << Date() << ": Loading reads from: " << READS_EC[i] << endl;
    reads.ReadAll( run_dir + "/" +  READS_EC[i] + ".fastb", True );
  }

  // determine the median of the read lengths
  int median_rl = -1;
  {
    vec <int> read_lengths( reads.size(), 0);
    for ( size_t i = 0; i < reads.size(); i++ ) 
      read_lengths[i] = reads[i].size();
    Sort( read_lengths );
    median_rl = Median( read_lengths );
  }

  VecPdfEntryVec cn_pdfs;
  vec<int64_t> n_kmer_hits;
  vec<double> CN_raw( unibases.size() );
  //
  // find the largest supported K so that K < median_rl
  // !! Please note that if K2 is significantly less than K, the calculated n_kmer_hits will 
  // not be trustable. Copy number estimations may also be affected.
  //
  #define ASSIGN_K( K, valname )  \
     if( K < median_rl ) valname = K;  
  int val = 0 ;
  FOR_ALL_K( ASSIGN_K,   val);
  int K2 = Min( val, K );
  ForceAssertGt( K2, 0 );
  if ( K2 < K ) {
    cout << Date() << ": *** Use  K= " << K2 << " for copy number calculation " << endl;
  }
  UnibaseCopyNumber3Core( K2, reads, unibases, NUM_THREADS, merged_head, unibases_head, 
            PLOIDY, VERBOSE, cout, cn_pdfs, CN_raw, n_kmer_hits, THRESH, ERR_RATE, CORRECT_BIAS );
  
     // Start of code to lower copy numbers.  First, don't do it if PLOIDY > 1 or
     // if we haven't been asked to lower copy numbers.

     if ( PLOIDY > 1 ) LOWER = False;
     if ( !LOWER )
     {    if (WRITE)
          {    cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS 
	             + ".predicted_count.k" + ToString(K)).c_str() );
               BinaryWrite3( run_dir + "/" + READS + "." + UNIPATHS + ".kmer_hits.k" 
	 	     + ToString(K), n_kmer_hits );
               BinaryWrite3( run_dir + "/" + READS + "." + UNIPATHS + ".cn_raw.k" 
	 	     + ToString(K), CN_raw );    }
          cout << Date( ) << ": UnibaseCopyNumber3 Done!" << endl;
          return 0;    }


     // Now load unipath placements on reference if needed.

     VecPlacementVec placements;
     if ( LOWER_VERBOSITY >= 3 )
     {    placements.ReadAll( run_dir + "/" + READS + ".unipaths.k" + ToString(K)
	       + ".locs" );    }
  
     // Load innie stats.
  
     vec<int> innie_sep, innie_dev;
     vec<double> innie_percent;
     fast_ifstream iin( nonjumps_file );
     String line;
     while(1)
     {    getline( iin, line );
          if ( iin.fail( ) ) break;
          istrstream iline( line.c_str( ) );
          int sep, dev;
          double iper;
          iline >> sep >> dev >> iper;
          innie_sep.push_back(sep), innie_dev.push_back(dev);
          innie_percent.push_back(iper);
          if ( LOWER_VERBOSITY >= 1 )
          {    cout << Date( ) << ": jump innies " << sep << " +/- " << dev << " (" 
	            << setiosflags(ios::fixed) << setprecision(1)
	            << iper << "%)" << endl;    }    }
  
     // Align jumping reads.
     
     cout << Date( ) << ": " 
          << "-----------------------------------------------------" << endl;

     // Set up ancillary data structures for unibases.
  
     size_t nuni = unibases.size( );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc, K );

     vec<double> GC(nuni);
     vec<int> ulen(nuni), ids( nuni, vec<int>::IDENTITY );
     for ( size_t i = 0; i < nuni; i++ )
          ulen[i] = unibases[i].size( );
     ParallelReverseSortSync( ulen, ids );
     for ( size_t i = 0; i < nuni; i++ )
          ulen[i] = unibases[i].size( );
     for ( size_t u = 0; u < nuni; u++ )
     {    if ( (int) u > to_rc[u] ) continue;
          int gc = 0;
          for ( int j = 0; j < unibases[u].isize( ); j++ )
          if ( unibases[u][j] == 1 || unibases[u][j] == 2 ) gc++;
          GC[u] = GC[ to_rc[u] ] 
              = double(gc) / double( unibases[u].isize( ) );    }

     // Compute unibase copy numbers.  

     // The following module will do following:
     // 1.
     // For each pair of unipaths that are connected by two or more links, predict
     // their order and separation.  Both orders may be possible.  Note that a more 
     // wholistic approach may be needed.
     // 2. 
     // Make the list of putative unipath links into a set of condensed links.
     // Here we use GapStats to merge links, and we apply a min_links_initial 
     // threshold.  A condensed links is accepted if its predicted overlap is no 
     // more than K-1, plus slop.
     vec<ulink_with_uids> condensed_links;
     vec<int> most_links_from( nuni, 0 ), most_links_to( nuni, 0 );
     vec< vec< triple<ho_interval,int,int> > > cov_from(nuni), cov_to(nuni);
     vec<double> cn_from(nuni, 0), cn_to(nuni, 0);

     cout << Date( ) << ": loading jumps" << endl;
     String temp_file = run_dir + "/tmp/"+ READS + ".temp_jumps.fastb";
     {    vecbasevector jreads( jreads_file );
	  for ( size_t i = 0; i < jreads.size( ); i++ )
	       jreads[i].ReverseComplement( );
	  cout << Date( ) << ": total jump reads = " 
	       << ToStringAddCommas( jreads.size( ) ) << endl;
	  for ( size_t i = 0; i < jreads.size( ); i++ )
	       jreads[i].resize(20);
	  jreads.WriteAll(temp_file);    
     }
     vec< triple<int64_t,int64_t,int> > JALIGNS;
     {    vec< triple<int64_t,int64_t,int> > jaligns;
	  SearchFastb2( temp_file, uni_file, 20, &jaligns, 0, MAX_PLACEMENTS );
	  for ( size_t i = 0; i < jaligns.size( ); i++ )
	       if ( jaligns[i].third >= 0 ) JALIGNS.push_back( jaligns[i] );    
     }
     // Directly convert alignments into segments, then get rid of alignments.
     vec<segalign> JSEGS;
     JSEGS.resize( JALIGNS.size( ) );
     for ( size_t i = 0; i < JALIGNS.size( ); i++ ) 
     {    JSEGS[i] = segalign( True, JALIGNS[i].first, 0, JALIGNS[i].second,
	       JALIGNS[i].third );    }
     Destroy(JALIGNS);
     ParallelSort(JSEGS, cmp_ridx);

     GetUnibaseSeps( run_dir, JUMP_READS, K, to_rc, JSEGS, unibases, innie_sep, 
	       innie_dev, innie_percent, CN_raw, // input 
	       condensed_links, cn_from, cn_to, cov_from, cov_to, most_links_from, 
	       most_links_to );

     // Re-evaluate the gaps using gap remodeling tools
     // Only gaps in condensed_links will be re-evaluated if REGAP is true
     vec <ulink_with_uids> condensed_links_save = condensed_links;
     if (  EXP_GAP_REMODEL ) {
          UnibaseGapRemodeling ( run_dir, JUMP_READS, unibases, K, to_rc, CN_raw, //input
		    condensed_links, cn_from, cn_to, cov_from, cov_to, most_links_from, 
		    most_links_to, //output
		    True ,  // REGAP mode ?
		    LOWER_VERBOSITY );
     }
     CheckGaps(to_rc,  condensed_links ); // check if gap(u1,u2) is the same as gap(u2*, u1*)

     map< pair<int,int>, int > record_gaps;
     if ( LOWER_VERBOSITY >= 2 )
     {    const double max_offby = 50.0;
          int total_links = 0, bad_links = 0;
          int nuni = unibases.size( );
          vec< vec<ulink_with_uids> > clinks(nuni);
          for ( size_t i = 0; i < condensed_links.size( ); i++ )
               clinks[ condensed_links[i].u1 ].push_back( condensed_links[i] );
          for ( int u1 = 0; u1 < nuni; u1++ )
          for ( int i = 0; i < clinks[u1].isize( ); i++ )
          {    const ulink_with_uids& l1 = clinks[u1][i];
               int u2 = l1.u2;
               //if ( !( u1 < u2 ) ) continue;
               Bool good = False;
               cout << "\n";
               for ( int pass = 1; pass <= 2; pass++ )
               {    ulink_with_uids l;
                    if ( pass == 1 ) l = l1;
                    else
                    {    Bool found = False;
                         for ( int j = 0; j < clinks[u2].isize( ); j++ )
                         {    if ( clinks[u2][j].u2 == u1 )
                              {    l = clinks[u2][j];
                                   found = True;
                                   break;    }    }    
                         if ( !found ) continue;    }
                    int u1 = l.u1, u2 = l.u2;
                    cout << u1 << "[CN_raw=" << CN_raw[u1] << ",len=" 
                         << unibases[u1].size( ) << "] --> " << u2 << "[CN_raw=" 
                         << CN_raw[u2] << ",len=" << unibases[u2].size( ) << "]: ";
                    cout << l.sep << " +/- " << l.dev << ", start1 = " << l.start1 
                         << ", stop2 = " << l.stop2 << ", nlinks = " << l.nlinks;
                    if ( LOWER_VERBOSITY >= 3 )
                    {    if ( placements[u1].empty( ) || placements[u2].empty( ) )
                         {    cout << ", not both mapped\n";
                              continue;    }
                         double infinity = 1000000000.0;
                         int best_sep = 0;
                         double best_offby = infinity;
                         for ( size_t j1 = 0; j1 < placements[u1].size(); j1++ )
                         {    for ( size_t j2 = 0; j2 < placements[u2].size(); j2++ )
                              {    const placement& p1 = placements[u1][j1];
                                   const placement& p2 = placements[u2][j2];
                                   if ( p1.Orient( ) != p2.Orient( ) ) continue;
                                   if ( p1.GenomeId( ) != p2.GenomeId( ) ) continue;
                                   int sep;
                                   double offby;
                                   if ( p1.Fw( ) ) sep = p2.pos( ) - p1.Pos( );
                                   else sep = p1.pos( ) - p2.Pos( );
                                   offby = double( Abs(l.sep-sep) ) / double(l.dev);
                                   if ( offby < best_offby )
                                   {    best_sep = sep;
                                        best_offby = offby;    }    }    }
                         if ( best_offby < max_offby )
                         {    cout << ", actual sep = " << best_sep << " (off by " 
                                   << best_offby << " devs)";    
                              record_gaps[ make_pair(u1,u2) ] = best_sep;
                              record_gaps[ make_pair(to_rc[u2],to_rc[u1]) ] = best_sep;
                              good = True;    }    }
                    cout << endl;    }
               if ( LOWER_VERBOSITY >= 3 )
               {    if ( placements[u1].empty( ) || placements[u2].empty( ) ) 
                         continue;
                    total_links++;
                    if ( !good ) bad_links++;    }    }
          if ( LOWER_VERBOSITY >= 3 )
          {    cout << "\nbad links = " << bad_links << ", total links = "
                    << total_links << ", "
                    << PERCENT_RATIO( 2, bad_links, total_links )
                    << " of links are bad" << endl;    }    }


     // Output the gap size estimations
     if ( WRITE_GAP ) {
          String outfile = run_dir + "/" + READS + "." + UNIPATHS 
               + ( WRITE_GAP_SUFFIX == "" ? "" : "." + WRITE_GAP_SUFFIX )
               + ".predicted_gap.k" + ToString(K);
          cout << Date( ) << ": Output gap size estimation to " + outfile << endl;
          ofstream out( outfile.c_str() );
          out << "#u1 u2 sep dev nlinks " << endl;
          for( int i =0; i < condensed_links.isize(); i++ ){
               const ulink_with_uids& link = condensed_links[i];
	       map< pair<int,int>, int > old_gaps;
	       for( int i =0; i < condensed_links.isize(); i++ ){
		    const ulink_with_uids& link = condensed_links_save[i]; 
		    old_gaps[ make_pair( link.u1, link.u2 ) ] = link.sep;
	       }
	       if ( link.sep == old_gaps[make_pair(link.u1,link.u2)] && GAP_STATS_REMODEL) continue;
               out << link.u1 << " " << link.u2 << " " << link.sep << " " << link.dev << " " << link.nlinks << endl;
          }
          out.close();
     }

     const double cn_thresh = 0.7;
     const int min_links_mult = 3;
     for ( size_t u = 0; u < nuni; u++ )
     {    int cnx;
          GetMostLikelyValue( cnx, cn_pdfs[u] );
          if ( cnx == 1 ) continue;
          if ( cn_from[u] >= cn_thresh * CN_raw[u]
	       && cn_to[u] >= cn_thresh * CN_raw[u] )
          {    int max_from = 0, max_to = 0;
               for ( int j = 0; j < cov_from[u].isize( ); j++ )
	       max_from = Max( max_from, cov_from[u][j].second );
               for ( int j = 0; j < cov_to[u].isize( ); j++ )
	            max_to = Max( max_to, cov_to[u][j].second );
               Bool overlap = False;
               for ( int j1 = 0; j1 < cov_from[u].isize( ); j1++ )
               {    if ( cov_from[u][j1].second < max_from/min_links_mult ) continue;
	            for ( int j2 = j1 + 1; j2 < cov_from[u].isize( ); j2++ )
                    {    if ( cov_from[u][j2].second < max_from/min_links_mult ) 
	                      continue;
	                 if ( Meets( cov_from[u][j1].first, cov_from[u][j2].first ) )
	                      overlap = True;    }    }
               for ( int j1 = 0; j1 < cov_to[u].isize( ); j1++ )
               {    if ( cov_to[u][j1].second < max_to/min_links_mult ) continue;
	            for ( int j2 = j1 + 1; j2 < cov_to[u].isize( ); j2++ )
                    {    if ( cov_to[u][j2].second < max_to/min_links_mult ) 
	                      continue;
	                 if ( Meets( cov_to[u][j1].first, cov_to[u][j2].first ) )
	                      overlap = True;    }    }
               if (overlap) continue;
               if ( LOWER_VERBOSITY >= 1 && (int) u <= to_rc[u] ) 
               {    cout << "\nsetting copy number of unipath " << u << "[l=" 
                         << unibases[u].isize( ) - K + 1 << ", gc=" 
                         << setiosflags(ios::fixed) << setprecision(1) 
                         << 100.0 * GC[u] << "%] to 1\n";
	            PRINT3( CN_raw[u], cn_from[u], cn_to[u] );    
	            for ( int j = 0; j < cov_to[u].isize( ); j++ )
                    {    cout << cov_to[u][j].first.Start( ) << " - "
                              << cov_to[u][j].first.Stop( ) << " [n=" 
	                      << cov_to[u][j].second << "]"
                              << " (u=" << cov_to[u][j].third << ")\n";    }    
	            for ( int j = 0; j < cov_from[u].isize( ); j++ )
                    {    cout << cov_from[u][j].first.Start( )
                              << " - " << cov_from[u][j].first.Stop( ) << " [n=" 
	                      << cov_from[u][j].second << "]"
                              << " (u=" << cov_from[u][j].third << ")\n";    }    }
               PdfEntryVec copyno_prob;
               copyno_prob.push_back( pdf_entry( 1, 1 ) );
               cn_pdfs[u] = copyno_prob;    }    }


     // gap estimation accuracy evaluation table
     if ( LOWER_VERBOSITY >= 3 ) {
          bool StatRemodelOnly = GAP_STATS_REMODEL;
	  map< pair<int,int>, int > old_gaps;
          for( int i =0; i < condensed_links.isize(); i++ ){
	    const ulink_with_uids& link = condensed_links_save[i]; 
	    old_gaps[ make_pair( link.u1, link.u2 ) ] = link.sep;
	  }
          cout << "Found " << record_gaps.size() << " gaps in reference genome. " << endl;
          cout << "Gap statistics " << endl;
	  if ( StatRemodelOnly )
	    cout << "Check only remodeled gaps " << endl;
          cout << "-------------------------------------------------------------" << endl;
          int off1 = 0, off2 = 0, off3 = 0, off4 = 0, off5 = 0, off5plus = 0;
          int bad_links = 0;
          int n = condensed_links.size();

          for( int i =0; i < condensed_links.isize(); i++ ){
               const ulink_with_uids& link = condensed_links[i];
	       if ( link.sep == old_gaps[make_pair(link.u1,link.u2)] && StatRemodelOnly ) {
		    n--;
		    continue;
	       }
               if ( record_gaps.find( make_pair(link.u1, link.u2) ) == record_gaps.end() ) {
                    bad_links++;
                    continue;
               }
               int real_gap = record_gaps[ make_pair(link.u1, link.u2) ] ;
               double dev = fabs( double( link.sep - real_gap ) /double( link.dev ) );
               if      ( dev < 1 ) off1++;
               else if ( dev < 2 ) off2++;
               else if ( dev < 3 ) off3++;
               else if ( dev < 4 ) off4++;
               else if ( dev < 5 ) off5++;
               else                off5plus++;
          }

          vec< vec<String> > table;
          table.push_back( MkVec<String>(  "dev", "#gaps", "%" ) );
          table.push_back( MkVec<String>(  "0 - 1", ToString(off1), ToString( off1 * 100/n) ) );
          table.push_back( MkVec<String>(  "1 - 2", ToString(off2), ToString( off2 * 100/n) ) );
          table.push_back( MkVec<String>(  "2 - 3", ToString(off3), ToString( off3 * 100/n) ) );
          table.push_back( MkVec<String>(  "3 - 4", ToString(off4), ToString( off4 * 100/n) ) );
          table.push_back( MkVec<String>(  "4 - 5", ToString(off5), ToString( off5 * 100/n) ) );
          table.push_back( MkVec<String>(  ">=  5", ToString(off5plus), ToString( off5plus * 100/n) ) );
          table.push_back( MkVec<String>("wrong gap", ToString(bad_links), ToString( bad_links * 100/n) ) );
          table.push_back( MkVec<String>("Total",ToString(n), "100" ) );
          PrintTabular(cout, table, 2, "rrr" );
     }


     if (WRITE)
     {    cn_pdfs.WriteAll( (run_dir + "/" + READS + "." + UNIPATHS 
               + ".predicted_count.k" + ToString(K)).c_str() );
          BinaryWrite3( run_dir + "/" + READS + "." + UNIPATHS + ".kmer_hits.k" 
	       + ToString(K), n_kmer_hits );
          BinaryWrite3( run_dir + "/" + READS + "." + UNIPATHS + ".cn_raw.k" 
               + ToString(K), CN_raw );    }
     cout << Date( ) << ": UnibaseCopyNumber3 Done!" << endl;
     return 0;    }
