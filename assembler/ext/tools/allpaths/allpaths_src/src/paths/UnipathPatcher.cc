///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// UnipathPatcher.  Place reads on the unipaths then attempt to close gaps.
// This shares code with UnipathFixer.   Compare with CloseUnipathGaps.
//
// Don't feed this EcoP15I reads without changing the code appropriately.
//
// Note that this code does not generate bona fide unipaths as output unless
// you ask it to.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

// MakeDepend: dependency PatcherCottage
// MakeDepend: dependency QueryLookupTable

#include <omp.h>
#include <sys/mman.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "FindGaps.h"
#include "Intvector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "ParseSet.h"
#include "Qualvector.h"
#include "VecTemplate.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "paths/GetNexts.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseUtils.h"
#include "paths/Ulink.h"
#include "paths/Unipath.h"
#include "paths/UnipathFixerTools.h"
#include "system/SharedMem.h"

void GapStatsAlt( vec<int> gap, vec<int> gapdev, int& gap_ave, int& gapdev_ave )
{
     // If there are less than six gaps, we directly compute their mean.
     // Otherwise, we attempt to remove outliers, as follows.  We sort 
     // the gaps and extract the middle half.  From this middle half, we compute the 
     // mean and standard deviation.  Then we select those gaps lying withing 5 
     // standard deviations of the mean, and form their mean.  

     vec<normal_distribution> S;
     if ( gap.size( ) >= 6 )
     {    vec<int> mid_gaps;
          SortSync( gap, gapdev );
          for ( unsigned int i = gap.size( )/4; i < 3*(1+gap.size( ))/4; i++ )
               mid_gaps.push_back( gap[i] );
          float sum1 = 0, sum2 = 0;
          for ( unsigned int i = 0; i < mid_gaps.size( ); i++ )
          {    sum1 += mid_gaps[i];
               sum2 += float(mid_gaps[i]) * float(mid_gaps[i]);    }
          float n = mid_gaps.size( );
          float mean = sum1/n;
          float sd = sqrt(sum2/n - mean * mean);
          float start = mean - 5 * sd, stop = mean + 5 * sd;
          for ( unsigned int i = 0; i < gap.size( ); i++ )
          {    if ( start <= gap[i] && gap[i] <= stop )
                    S.push( gap[i], gapdev[i] );    }    }
     else
     {    for ( int l = 0; l < gap.isize( ); l++ )
               S.push( gap[l], gapdev[l] );    }
     normal_distribution s = CombineNormalDistributions(S);
     gap_ave = int(round(s.mu_));
     gapdev_ave = int(round(s.sigma_));    }

int main( int argc, char *argv[] )
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
          "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(IN_HEAD, "all_reads");
     CommandArgument_String_OrDefault(OUT_HEAD, "all_reads.patched");
     CommandArgument_String_OrDefault(JUMP_READS, "jump_reads_filt_cpd");
     CommandArgument_String_OrDefault(FRAG_READS, "frag_reads_filt_cpd");
     CommandArgument_String_OrDefault(FRAG_READS_EDIT, "frag_reads_edit");
     CommandArgument_Bool_OrDefault(PRINT_ALIGNMENTS, False);
     CommandArgument_Bool_OrDefault(PRINT_SEGMENTS, False);
     CommandArgument_Bool_OrDefault(SHOW_ALL, False);
     CommandArgument_Int_OrDefault(MAX_PLACEMENTS, 50);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_String_OrDefault_Doc(LOG, "",
          "{option_1,...,option_n} where option_i is one of\n"
          "ALL: all options\n"
          "SOME: many options\n"
          "LINKS: print links\n"
          "LINKS_DETAILED: print detailed links\n"
          "POSSIBLE_JOINS: print possible joins\n"
          "ATTEMPTED_JOINS: print attempted joins\n"
          "ASSEMBLY: log minimal assembly details\n"
          "ASSEMBLYi: log assembly details, level i, i in {1,2,3}\n"
          "ALIGN_READS: alignments of reads used\n"
          "ALIGNS: alignments for successful joins\n"
          "ALIGNS_ALL: alignments for attempted joins\n"
          "ACTION: log calls and returns during joining\n"
          "CORRECT: details of error correction during gap closing\n"
          "INTERESTING: generate dot file for links on interesting part of graph");
     CommandArgument_String_OrDefault_Doc(LR, "", 
          "a list of pairs \"(u1 u2)\" of unipaths to attempt to join "
          "(rather than pairs that are identified for joining)");
     CommandArgument_Bool_OrDefault(CHECKPOINT, False);
     CommandArgument_String_OrDefault_Doc(CHECKPOINT_HEAD, "",
          "head for granular checkpoints of alignment results");
     CommandArgument_String_OrDefault(PATCHDIR, "unipath_patch");
     CommandArgument_Bool_OrDefault(BONA_FIDE_UNIPATHS, False);
     CommandArgument_String_OrDefault_Doc(JOINS, "",
          "If specified, do only these joins.");
     CommandArgument_LongLong_OrDefault(MAX_MEMORY_GB,0);
     EndCommandArguments;

     // Thread control
   
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     SetMaxMemory(MAX_MEMORY_GB<<30);

     // Heuristic constants.

     const int min_links = 5;
     const int MAX_PATHS = 100;
     const int MAX_PATH_ITERATIONS = 100000;
     const int MAX_EDGES = 2000;
     const int MAX_READS = 1000;
     const int MAX_JOINS = 10;
     const int MIN_OVERLAP_END = 1000000000;

     // Start.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA, run_dir = PRE + "/" + DATA + "/" + RUN;
     cout << Date( ) << ": " << run_dir << endl;
     Mkdir777( run_dir + "/" + PATCHDIR );
     String temp_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.fastb";
     String unifile = run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K);
     if ( !IsRegularFile( unifile ) )
       FatalErr( "Can't find input file " << unifile );

     // Define logging options.

     vec<String> log_options, allowed_log_options;
     ParseStringSet( LOG, log_options );
     allowed_log_options.push_back( "ALL", "SOME", "LINKS", "LINKS_DETAILED" );
     allowed_log_options.push_back( "INTERESTING", "ACTION" );
     allowed_log_options.push_back( "POSSIBLE_JOINS", "ATTEMPTED_JOINS" );
     allowed_log_options.push_back( "ASSEMBLY", "CORRECT" );
     allowed_log_options.push_back( "ASSEMBLY1", "ASSEMBLY2", "ASSEMBLY3" );
     allowed_log_options.push_back( "ALIGN_READS", "ALIGNS", "ALIGNS_ALL" );
     Sort(allowed_log_options);
     if ( !BinSubset( log_options, allowed_log_options ) )
     {    cout << "Illegal log option." << endl << "Abort." << endl;
          exit(1);    }
     #define REQUESTED(x) Member( log_options, String(x) )
     Bool log_all = REQUESTED( "ALL" );
     Bool log_some = log_all || REQUESTED( "SOME" );
     Bool log_links = log_all || REQUESTED( "LINKS" );
     Bool log_links_detailed = log_all || REQUESTED( "LINKS_DETAILED" );
     Bool log_possible_joins = log_some || REQUESTED( "POSSIBLE_JOINS" );
     Bool log_attempted_joins = log_some || REQUESTED( "ATTEMPTED_JOINS" );
     Bool log_aligns = log_some || REQUESTED( "ALIGNS" );
     Bool log_action = log_all || REQUESTED( "ACTION" );
     Bool log_interesting = log_all || REQUESTED( "INTERESTING" );

     // Parse arguments.

     vec< pair<int,int> > LR_to_process;
     {    vec<String> LRs;
          ParseStringSet( LR, LRs );
          for ( int i = 0; i < LRs.isize( ); i++ )
          {    int L = LRs[i].Between( "(", " " ).Int( );
               int R = LRs[i].Between( " ", ")" ).Int( );
               LR_to_process.push( L, R );    }
          Sort(LR_to_process);    }
     vec<int> joins;
     ParseIntSet( JOINS, joins );

     // Load innie stats.

     vec<int> innie_sep, innie_dev;
     vec<double> innie_percent;
     fast_ifstream iin( run_dir + "/" + JUMP_READS + ".outies" );
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
          cout << Date( ) << ": jump innies " << sep << " +/- " << dev << " (" 
               << setiosflags(ios::fixed) << setprecision(1)
               << iper << "%)" << endl;    }

     // Align reads.

     vec< triple<int64_t,int64_t,int> > ALIGNS, JALIGNS;
     String ALIGNS_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.ALIGNS";
     String JALIGNS_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.JALIGNS";
     if ( CHECKPOINT && IsRegularFile(ALIGNS_file) && IsRegularFile(JALIGNS_file) )
     {    BinaryRead3( ALIGNS_file, ALIGNS );
          BinaryRead3( JALIGNS_file, JALIGNS );    }
     else
     {    double align_clock = WallClockTime( );
          AlignReadsToUnipaths( run_dir, JUMP_READS, FRAG_READS, FRAG_READS_EDIT, 
               True, MAX_PLACEMENTS, unifile, ALIGNS, JALIGNS, CHECKPOINT_HEAD );
          cout << Date( ) << ": " << TimeSince(align_clock)
               << " used aligning reads to unipaths" << endl;
          BinaryWrite3( ALIGNS_file, ALIGNS );
          BinaryWrite3( JALIGNS_file, JALIGNS );
          if (PRINT_ALIGNMENTS)
          {    for ( size_t i = 0; i < ALIGNS.size( ); i++ )
               {    cout << "fragment read " << ALIGNS[i].first << " --> "
                         << "unipath " << ALIGNS[i].second << " starting at " 
                         << ALIGNS[i].third << "\n";    }
               for ( size_t i = 0; i < JALIGNS.size( ); i++ )
               {    cout << "jumping read " << JALIGNS[i].first << " --> "
                         << "unipath " << JALIGNS[i].second << " starting at " 
                         << JALIGNS[i].third << "\n";    }    }    }

     // Load unibases.

     cout << Date( ) << ": load unibases" << endl;
     vecbasevector unibases( run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K) );
     size_t nuni = unibases.size( );
     for ( int i = 0; i < LR_to_process.isize( ); i++ )
     {    ForceAssertGe( LR_to_process[i].first, 0 );
          ForceAssertGe( LR_to_process[i].second, 0 );
          ForceAssertLt( LR_to_process[i].first, (int) nuni );
          ForceAssertLt( LR_to_process[i].second, (int) nuni );    }

     // Set up ancillary data structures for unibases.

     cout << Date( ) << ": computing nexts" << endl;
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     cout << Date( ) << ": computing unibase involution" << endl;
     vec<int> to_rc;
     String TORC_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.TORC";
     if ( CHECKPOINT && IsRegularFile(TORC_file) ) BinaryRead3( TORC_file, to_rc );
     else
     {    UnibaseInvolution( unibases, to_rc, K );
          if (CHECKPOINT) BinaryWrite3( TORC_file, to_rc );    }

     // Directly convert alignments into segments, then get rid of alignments.

     cout << Date( ) << ": converting alignments into segments" << endl;
     vec<segalign> SEGS, JSEGS;
     String SEGS_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.SEGS";
     String JSEGS_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.JSEGS";
     if ( CHECKPOINT && IsRegularFile(SEGS_file) ) 
     {    BinaryRead3( SEGS_file, SEGS );
          Destroy(ALIGNS);
          BinaryRead3( JSEGS_file, JSEGS );
          Destroy(JALIGNS);    }
     else
     {    SEGS.resize( ALIGNS.size( ) ), JSEGS.resize( JALIGNS.size( ) );
          for ( size_t i = 0; i < ALIGNS.size( ); i++ )
          {    SEGS[i] = segalign( True, ALIGNS[i].first, 0, ALIGNS[i].second,
                    ALIGNS[i].third );    }
          Destroy(ALIGNS);
          ParallelSort(SEGS);
          for ( size_t i = 0; i < JALIGNS.size( ); i++ )
          {    JSEGS[i] = segalign( True, JALIGNS[i].first, 0, JALIGNS[i].second,
                    JALIGNS[i].third );    }
          Destroy(JALIGNS);
          ParallelSort(JSEGS);
          BinaryWrite3( SEGS_file, SEGS );
          BinaryWrite3( JSEGS_file, JSEGS );    }

     // Index the aligned segments.

     cout << Date( ) << ": indexing segments" << endl;
     vec<size_t> S_START(nuni+1);
     {    size_t SEG_POS = 0;
          for ( size_t u = 0; u <= nuni; u++ )
          {    while( SEG_POS < SEGS.size( ) && SEGS[SEG_POS].u < u ) ++SEG_POS;
               S_START[u] = SEG_POS;    }    }
     vec<size_t> JS_START(nuni+1);
     {    size_t SEG_POS = 0;
          for ( size_t u = 0; u <= nuni; u++ )
          {    while( SEG_POS < JSEGS.size( ) && JSEGS[SEG_POS].u < u ) 
                    ++SEG_POS;
               JS_START[u] = SEG_POS;    }    }

     // Compute unibase copy numbers.  Derived with as few changes as possible from
     // UnibaseCopyNumber2.

     int PLOIDY = StringOfFile( run_dir + "/ploidy", 1 ).Int( );
     cout << Date( ) << ": computing unibase copy numbers, "
          << "ploidy = " << PLOIDY << endl;
     vec<int> CN;
     String CN_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.CN";
     if ( CHECKPOINT && IsRegularFile(CN_file) ) BinaryRead3( CN_file, CN );
     else
     {    
          // Quick and dirty calculation.

          CN.resize(nuni);
          vec<int> hits( nuni, 0 ), ulen(nuni), ids( nuni, vec<int>::IDENTITY );
          for ( size_t i = 0; i < SEGS.size( ); i++ )
          {    hits[ SEGS[i].u ]++;
               hits[ to_rc[ SEGS[i].u ] ]++;    }
          for ( size_t i = 0; i < nuni; i++ )
               ulen[i] = unibases[i].size( );
          ParallelReverseSortSync( ulen, ids );
          int sample = Min( 100, (int) nuni );
          int64_t total_bases = 0, total_hits = 0;
          for ( int i = 0; i < sample; i++ )
          {    total_bases += ulen[i]; total_hits += hits[ ids[i] ];    }
          for ( size_t i = 0; i < nuni; i++ )
               ulen[i] = unibases[i].size( );
          double hits_per_base = double(total_hits) / double(total_bases);
          for ( size_t u = 0; u < nuni; u++ )
          {    if ( (int) u > to_rc[u] ) continue;
               double cov = ( double(hits[u]) / double(ulen[u]) ) / hits_per_base;
               CN[u] = int( round( double(PLOIDY) * cov ) );
               CN[ to_rc[u] ] = int( round( double(PLOIDY) * cov ) );
               /*
               cout << "u = " << u << ", ulen = " << ulen[u]
                    << ", cov = " << cov << ", cn = " << CN[u];
               if ( CN[u] != 1 ) cout << " ***********************************";
               cout << "\n";    
               */
               }

          // No clue why, but memory usage blew up (should retry):

          /*
          UnibaseCopyNumbersFromSegAligns( K, unibases, to_rc, PLOIDY, run_dir, 
               FRAG_READS, SEGS, S_START, CN );
          */

          if (CHECKPOINT) BinaryWrite3( CN_file, CN );    }

     // Heuristics.

     const double max_devs = 6.0;
     const int min_kmers = 120;
     const int trans_depth = 3;
     const int max_sep_for_trans_del = 200;
     const int min_links_initial = 2;
     const int min_glue = 300;
     const int bi_mult = 5;
     const double min_innie_percent = 5.0;

     // For each pair of unipaths that are connected by two or more links, predict
     // their order and separation.  Both orders may be possible.  Note that a more 
     // wholistic approach may be needed.

     vec<ulink_with_uids> ulinks;
     MakeUlinks( run_dir, PATCHDIR, FRAG_READS, JUMP_READS, K, PLOIDY, CN, to_rc, 
          SEGS, JSEGS, unibases, innie_sep, innie_dev, innie_percent, 
          min_innie_percent, log_links_detailed, min_kmers, max_devs, ulinks, 
          CHECKPOINT );

     // Make the list of putative unipath links into a set of condensed links.
     // Here we use GapStats to merge links, and we apply a min_links_initial 
     // threshold.  A condensed links is accepted if its predicted overlap is no 
     // more than K-1, plus slop.

     cout << Date( ) << ": building condensed links" << endl;
     vec<ulink_with_uids> condensed_links;
     for ( size_t i = 0; i < ulinks.size( ); i++ )
     {    size_t j;
          int u1 = ulinks[i].u1, u2 = ulinks[i].u2;
          ForceAssertLe( CN[u1], PLOIDY ); 
          ForceAssertLe( CN[u2], PLOIDY );
          for ( j = i + 1; j < ulinks.size( ); j++ )
               if ( ulinks[j].u1 != u1 || ulinks[j].u2 != u2 ) break;
          if ( (int) (j - i) >= min_links_initial )
          {    vec<int> seps, devs;
               int sep, dev, start1 = unibases[u1].size( ), stop2 = 0;
               for ( size_t m = i; m < j; m++ )
               {    const ulink_with_uids& f = ulinks[m];
                    seps.push_back(f.sep), devs.push_back(f.dev);    
                    start1 = Min( start1, f.start1 );
                    stop2 = Max( stop2, f.stop2 );    }
               GapStatsAlt( seps, devs, sep, dev );
               if ( sep + max_devs * dev >= -(K-1) )
                    condensed_links.push( u1, u2, sep, dev, start1, stop2, j-i );   }
          i = j - 1;    }
     if (log_links)
     {    for( size_t i = 0; i < condensed_links.size( ); i++ )
          {    const ulink_with_uids& l = condensed_links[i];
               cout << "\n" << l.u1 << "[CN=" << CN[l.u1] << ",len=" 
                    << unibases[l.u1].size( ) << "] --> " << l.u2 << "[CN=" 
                    << CN[l.u2] << ",len=" << unibases[l.u2].size( ) << "]: ";
               cout << l.sep << " +/- " << l.dev << ", start1 = " << l.start1 
                    << ", stop2 = " << l.stop2 << ", nlinks = " << l.nlinks
                    << endl;    }    }
     for( size_t i = 0; i < condensed_links.size( ); i++ )
     {    const ulink_with_uids& l = condensed_links[i];
          ForceAssertLe( CN[l.u1], 2 );
          ForceAssertLe( CN[l.u2], 2 );    }

     // Build a digraph out of the condensed links.

     cout << Date( ) << ": building digraph from condensed links" << endl;
     digraphE<ulink> G;
     {    vec< vec<int> > from(nuni), to(nuni);
          vec< vec<int> > from_edge_obj(nuni), to_edge_obj(nuni);
          vec<ulink> edges;
          for( size_t i = 0; i < condensed_links.size( ); i++ )
          {    const ulink_with_uids& l = condensed_links[i];
               int u1 = l.u1, u2 = l.u2;
               from[u1].push_back(u2);
               from_edge_obj[u1].push_back( edges.size( ) );
               to[u2].push_back(u1);
               to_edge_obj[u2].push_back( edges.size( ) );
               edges.push( l.sep, l.dev, l.start1, l.stop2, l.nlinks );    }
          for ( size_t u = 0; u < nuni; u++ )
          {    UniqueSortSync( from[u], from_edge_obj[u] );
               UniqueSortSync( to[u], to_edge_obj[u] );    }
          G.Initialize( from, to, edges, to_edge_obj, from_edge_obj );    }

     // Reduce memory usage.

     Destroy(ulinks), Destroy(condensed_links), Destroy(CN);

     // Start cleaning up graph.

     cout << Date( ) << ": cleaning up graph" << endl;
     DeleteBi( G, bi_mult, log_links );
     DeleteBiAnother( G, unibases, log_links );

     // Look for indirect links.

     cout << Date( ) << ": looking for indirect links 1" << endl;
     vec< triple<int,int,ulink> > new_links;
     for ( int c = 0; c < G.N( ); c++ )
     {    for ( int j1 = 0; j1 < G.To(c).isize( ); j1++ )
          {    int a = G.To(c)[j1];
               const ulink& ea = G.EdgeObjectByIndexTo( c, j1 );
               for ( int j2 = 0; j2 < G.To(c).isize( ); j2++ )
               {    int b = G.To(c)[j2];
                    if ( a == b ) continue;
                    const ulink& eb = G.EdgeObjectByIndexTo( c, j2 );
                    int sep_ab = ea.Sep( ) - eb.Sep( ) - unibases[b].isize( );
                    int dev_ab = int( round( sqrt( 
                         ea.Dev( )*ea.Dev( ) + eb.Dev( )*eb.Dev( ) ) ) );
                    if ( sep_ab + max_devs * dev_ab < -(K-1) ) continue;
                    ulink x( sep_ab, dev_ab, ea.start1, unibases[b].size( ),
                         Min( ea.nlinks, eb.nlinks ) );
                    if (log_links)
                    {    cout << "new link " << a << " --> " << b
                              << " (" << sep_ab << " +/- " << dev_ab << ")"
                              << " from " << a << " --> " << c << " and "
                              << b << " --> " << c << "\n";    }
                    new_links.push( a, b, x );    }    }    }
     cout << Date( ) << ": looking for indirect links 2" << endl;
     for ( int c = 0; c < G.N( ); c++ )
     {    for ( int j1 = 0; j1 < G.From(c).isize( ); j1++ )
          {    int a = G.From(c)[j1];
               const ulink& ea = G.EdgeObjectByIndexFrom( c, j1 );
               for ( int j2 = 0; j2 < G.From(c).isize( ); j2++ )
               {    int b = G.From(c)[j2];
                    if ( a == b ) continue;
                    const ulink& eb = G.EdgeObjectByIndexFrom( c, j2 );
                    int sep_ab = eb.Sep( ) - ea.Sep( ) - unibases[a].isize( );
                    int dev_ab = int( round( sqrt( 
                         ea.Dev( )*ea.Dev( ) + eb.Dev( )*eb.Dev( ) ) ) );
                    if ( sep_ab + max_devs * dev_ab < -(K-1) ) continue;
                    ulink x( sep_ab, dev_ab, ea.start1, unibases[b].size( ),
                         Min( ea.nlinks, eb.nlinks ) );
                    if (log_links)
                    {    cout << "new link " << a << " --> " << b
                              << " (" << sep_ab << " +/- " << dev_ab << ")"
                              << " from " << c << " --> " << a << " and "
                              << c << " --> " << b << "\n";    }
                    new_links.push( a, b, x );    }    }    }
     cout << Date( ) << ": looking for indirect links 3" << endl;
     for ( int i = 0; i < new_links.isize( ); i++ )
     {    int a = new_links[i].first;
          int b = new_links[i].second;
          const ulink& x = new_links[i].third;
          Bool done = False;
          for ( int l = 0; l < G.From(a).isize( ); l++ )
          {    int d = G.From(a)[l];
               if ( d == b )
               {    ulink& e = G.EdgeObjectByIndexFromMutable( a, l );
                    e = CombineUlinks( x, e );
                    done = True;
                    break;    }    }
          if ( !done ) G.AddEdge( a, b, x );    }

     // Clean up the graph.

     cout << Date( ) << ": calling DeleteNlinks" << endl;
     DeleteNlinks( G, min_links, log_links );
     cout << Date( ) << ": calling DeleteTrans" << endl;
     DeleteTrans( G, unibases, trans_depth, max_sep_for_trans_del, log_links );
     cout << Date( ) << ": calling DeleteWeak" << endl;
     DeleteWeak( G, unibases, min_glue, log_links );
     cout << Date( ) << ": removing dead edge objects" << endl;
     G.RemoveDeadEdgeObjects( );

     // Define the interesting part of the graph and dump dot file.  Note bad use
     // of temp files.

     if (log_interesting)
     {    cout << Date( ) << ": defining interesting part of graph" << endl;
          String temp_fileb = "xxxxx.fastb", temp_filea = "xxxxx.aligns";
          String QLT = "QueryLookupTable" + ARG(K, 12) + ARG(MM, 12) + ARG(MC, 0.15)
               + ARG(SMITH_WAT, True) + ARG(L, data_dir + "/../genome.lookup")
               + ARG(SEQS, temp_fileb) + ARG(PARSEABLE, True);
          unibases.WriteAll(temp_fileb);
          System( QLT + " > " + temp_filea );
          vec<look_align> aligns;
          LoadLookAligns( temp_filea, aligns );
          Remove(temp_fileb), Remove(temp_filea);
          vec<String> orientation( nuni, "unknown" );
          for ( int i = 0; i < aligns.isize( ); i++ )
          {    int fw = aligns[i].Fw1( );
               int u = aligns[i].query_id;
               if ( fw && orientation[u] == "unknown" ) orientation[u] = "fw";
               if ( fw && orientation[u] == "rc" ) orientation[u] = "mixed";
               if ( !fw && orientation[u] == "unknown" ) orientation[u] = "rc";
               if ( !fw && orientation[u] == "fw" ) orientation[u] = "mixed";    }
          vec< vec<int> > comp;
          G.Components(comp);
          vec<int> use;
          for ( int i = 0; i < comp.isize( ); i++ )
          {    int fw = 0, rc = 0, total = 0;
               for ( int j = 0; j < comp[i].isize( ); j++ )
               {    int u = comp[i][j];
                    if ( orientation[u] == "fw" ) fw += unibases[u].size( );
                    if ( orientation[u] == "rc" ) rc += unibases[u].size( );
                    total += unibases[u].size( );    }
               if ( double(rc)/double(total) >= 0.9 ) continue;
               if ( comp[i].solo() && unibases[ comp[i][0] ].isize( ) - K + 1 < 500 )
                    continue;
               use.append( comp[i] );    }
          Sort(use);
          Ofstream( out, "xxx.dot" );
          vec<String> labels( G.N( ) );
          for ( int i = 0; i < G.N( ); i++ )
               labels[i] = ToString(i);
          G.DOT_vl( out, labels, use );    }

     // Print edges.

     if (log_links)
     {    cout << "\nGRAPH EDGES:\n";
          for ( int u1 = 0; u1 < G.N( ); u1++ )
          {    for ( int j = 0; j < G.From(u1).isize( ); j++ )
               {    int u2 = G.From(u1)[j];
                    const ulink& e = G.EdgeObjectByIndexFrom( u1, j );
                    cout << u1 << "[len=" << unibases[u1].size( ) << "] --> " << u2 
                         << "[len=" << unibases[u2].size( ) << "]: "
                         << e.sep << " +/- " << e.dev << ", start1 = " << e.start1 
                         << ", stop2 = " << e.stop2 << ", join length = " 
                              << unibases[u1].isize( ) - e.start1 + e.stop2
                         << ", links = " << e.nlinks << "\n";    }    }
          cout << "\n";    }

     // For every edge, try to figure out if there is bridge across it.

     cout << Date( ) << ": looking for bridges between unipaths" << endl;
     double bridge_clock = WallClockTime( );
     const int max_bridge_passes = 100;
     const int max_bridge_pile = 10000;
     const int allowed_overlap = 200;
     vec<Bool> bridge( G.EdgeObjectCount( ) );
     int no_bridges = 0, overpile = 0, bridge_impossible = 0;
     for ( size_t u1 = 0; u1 < nuni; u1++ )
     {    for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               if (log_links)
                    cout << "\nlooking from " << u1 << " to " << u2 << endl;
               int e = G.EdgeObjectIndexByIndexFrom( u1, j );
               const ulink& E = G.EdgeObjectByIndexFrom( u1, j );
               int sep = E.sep, dev = E.dev;
               if ( sep < -K + 1 - 3*dev - allowed_overlap )
               {    if (log_links)
                    {    cout << "not possible, sep/dev = " << sep << " +/- "
                              << dev << endl;    }
                    bridge_impossible++;
                    continue;    }
               vec< pair<int,int> > theres; // (unipath, separation)
               theres.push( u1, 0 );
               int i = 0;
               for ( int pass = 1; pass <= max_bridge_passes; pass++ )
               {    int nt = theres.size( );
                    for ( ; i < nt; i++ )
                    {    int u = theres[i].first, usep = theres[i].second;
                         for ( int l = 0; l < nexts[u].isize( ); l++ )
                         {    int v = nexts[u][l];
                              if ( v == u2 )
                              {    bridge[e] = True;
                                   if (log_links) cout << "found bridge";
                                   goto tail;    }
                              int newsep = sep + unibases[v].isize( ) - K + 1;
                              if ( newsep > sep + 3*dev ) continue;
                              theres.push( v, newsep );    }    }
                    if ( theres.isize( ) == nt )
                    {    if (log_links) cout << "no bridge found" << endl;
                         no_bridges++;
                         break;    }
                    if ( theres.isize( ) >= max_bridge_pile )
                    {    if (log_links) cout << "bridge pile limit exceeded" << endl;
                         overpile++;
                         break;    }    }
               tail: continue;    }    }
     int nedges = G.EdgeObjectCount( ), bridged = Sum(bridge);
     int need_more_iterations = nedges - bridged - no_bridges
          - overpile - bridge_impossible;
     cout << Date( ) << ": bridged: " << bridged << ", unbridgable: "
          << no_bridges << ", too close: " << bridge_impossible << endl;
     cout << Date( ) << ": overpile: " << overpile
          << ", need more iterations: " << need_more_iterations << endl;
     cout << Date( ) << ": done bridging, time used = " 
          << TimeSince(bridge_clock) << endl;
     vec< pair<int,int> > to_process;
     vec< pair<int,int> > to_process_sepdev;
     for ( size_t u1 = 0; u1 < nuni; u1++ )
     {    for ( int j = 0; j < G.From(u1).isize( ); j++ )
          {    int u2 = G.From(u1)[j];
               int e = G.EdgeObjectIndexByIndexFrom( u1, j );
               if ( bridge[e] ) continue;
               const ulink& E = G.EdgeObjectByIndexFrom( u1, j );
               int sep = E.sep, dev = E.dev;
               if (log_links) 
               {    cout << "maybe process ";
                    PRINT4( u1, u2, sep, dev );        }
               if ( sep < 250 ) 
               {    to_process.push( u1, u2 );
                    to_process_sepdev.push( sep, dev );    }    }    }

     // Choose between pairs and their reverse complements.

     Sort(to_process);
     vec<Bool> to_delete( to_process.size( ), False );
     for ( int i = 0; i < to_process.isize( ); i++ )
     {    int u1 = to_process[i].first, u2 = to_process[i].second;
          int r1 = to_rc[u2], r2 = to_rc[u1];
          if ( make_pair(u1,u2) < make_pair(r1,r2)
               && BinMember( to_process, make_pair(r1,r2) ) )
          {    to_delete[i] = True;    }    }
     EraseIf( to_process, to_delete );
     EraseIf( to_process_sepdev, to_delete );

     // Fish for connections between dead ends.

     cout << Date( ) << ": looking for dead ends" << endl;
     vec<int> dead_fw, dead_rc;
     for ( int i = 0; i < to_process.isize( ); i++ )
     {    dead_fw.push_back( to_process[i].first );
          dead_rc.push_back( to_process[i].second );    }
     int nde = dead_fw.size( ) + dead_rc.size( );
     cout << Date( ) << ": found " << nde << " dead ends" << endl;

     // Sort stuff.

     cout << Date( ) << ": sorting SUGS" << endl;
     for ( int pass = 1; pass <= 2; pass++ )
     {    vec<segalign>& SUGS = ( pass == 1 ? SEGS : JSEGS );
          ParallelSort(SUGS);    }

     // Find linked dead ends.  Promiscuous links are introduced when links go over
     // stuff.  We should only use the closest links.

     vec< pair<int,int> > joiners;
     for ( int i = 0; i < to_process.isize( ); i++ )
          joiners.push( i, dead_fw.isize( ) + i );

     // Build join data.  Gather up the reads that are placed within 30 bases of a 
     // dead end, along with their partners.  Note that unpaired reads are ignored, 
     // which is not necessarily the right thing to do.  Also we always choose the
     // innie orientation for jumps, which is usually wrong.  

     String JOINDATA_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.JOINDATA";
     String JDLEN_file = run_dir + "/" + PATCHDIR + "/UnipathPatcher.JDLEN";
     vec<size_t> join_data_offsets;
     if ( CHECKPOINT && IsRegularFile(JOINDATA_file) && IsRegularFile(JDLEN_file) )
     {    Destroy(SEGS), Destroy(S_START), Destroy(JSEGS), Destroy(JS_START);
          Destroy(to_rc);
          BinaryReader::readFile(JDLEN_file.c_str(),&join_data_offsets);    }
     else
     {    BuildJoinData( dead_fw, dead_rc, to_process_sepdev, unibases, to_rc, 
               run_dir, PATCHDIR, FRAG_READS, JUMP_READS, SEGS, JSEGS, S_START, 
               JS_START, innie_sep, innie_dev, innie_percent, min_innie_percent,
               join_data_offsets, joiners, JOINDATA_file );
          if (CHECKPOINT)
              BinaryWriter::writeFile(JDLEN_file.c_str(),join_data_offsets);   }
     cout << Date( ) << ": join data created, memory usage = "
          << ToStringAddCommas( MemUsageBytes( ) ) << endl;

     // Try to make joins.

     vec< pair<int,int> > Xjoiners;
     for ( int i = 0; i < joiners.isize( ); i++ )
     {    int v1 = joiners[i].first, v2 = joiners[i].second;
          int u1 = dead_fw[v1], u2 = dead_rc[ v2 - dead_fw.isize( ) ];
          Xjoiners.push( u1, u2 );    }
     bvec3 new_stuff;
     vec< vec<int> > start_stop;
     vec<String> reports;
     double join_clock = WallClockTime( );
     int attempted_joins, npatches;
     vec<Bool> joined, perfect;
     MakeJoins( Xjoiners, unibases, "UnipathPatcher", K, data_dir, run_dir,
          PATCHDIR, attempted_joins, joined, npatches, perfect, new_stuff, 
          start_stop, reports, NUM_THREADS, joins, LOG, log_attempted_joins, log_action, 
          MAX_PATHS, MAX_PATH_ITERATIONS, MAX_EDGES, MAX_READS, MAX_JOINS, 
          MIN_OVERLAP_END, JOINDATA_file, join_data_offsets, LR_to_process );

     // Save results.

     cout << Date( ) << ": join complete" << endl;
     size_t nnew = 0;
     for ( size_t i = 0; i < new_stuff.size( ); i++ )
          nnew += new_stuff[i].size( );
     vecbasevector all_new;
     all_new.reserve(nnew);
     for ( size_t i = 0; i < new_stuff.size( ); i++ )
     {
         bvec3::const_reference vbv = new_stuff[i];
          all_new.append( vbv.begin(), vbv.end() );
     }
     cout << Date( ) << ": all_new constructed" << endl;
     if (CHECKPOINT) 
          all_new.WriteAll( run_dir + "/" + PATCHDIR + "/UnipathPatcher.NEW" );
     new_stuff.clear().shrink_to_fit();

     // Print reports.

     Bool have_report = False; 
     for ( int i = 0; i < reports.isize( ); i++ )
          if ( reports[i].size( ) > 0 ) have_report = True;    
     if (have_report) cout << "\n";
     for ( int i = 0; i < reports.isize( ); i++ )
          cout << reports[i];
     if (have_report) cout << "\n";
     cout << Date( ) << ": " << TimeSince(join_clock) << " used in joining" << endl;

     // Rebuild and write unipaths.

     if (WRITE)
     {    
          // Easy way.

          String KS = ToString(K);
          if ( !BONA_FIDE_UNIPATHS )
          {    unibases.Append(all_new);
               unibases.WriteAll( run_dir + "/" + OUT_HEAD + ".unibases.k" 
                    + KS );    }

          // Hard way.

          else
          {    // Rebuild unipaths.

               cout << Date( ) << ": rebuilding unipaths" << endl;
               double rebuild_clock = WallClockTime( );
               cout << Date( ) << ": " << unibases.size( ) << " old unibases" 
                    << endl;
               for ( size_t u = 0; u < nuni; u++ )
               {    for ( int j = 0; j < nexts[u].isize( ); j++ )
                    {    basevector b1 = unibases[u], b2 = unibases[ nexts[u][j] ];
                         b1.SetToSubOf( b1, b1.isize( ) - K, K );
                         b2.SetToSubOf( b2, K-1, 1 );
                         basevector b = Cat( b1, b2 );
                         unibases.push_back_reserve(b);    }    }
               size_t N = unibases.size( );
               unibases.Append(all_new);
               Destroy(all_new);
               cout << Date( ) << ": added " << unibases.size( ) - N
                    << " sequences" << endl;
               vecKmerPath paths, pathsrc, unipaths;
               vec<tagged_rpint> pathsdb, unipathsdb;
          
               // Run ReadsToPaths. On a 48-processor machine, using 32 processors
               // seems to be faster than using all 48, at least based on one data
               // point.  Therefore we use 2/3 of the threads.

               cout << Date( ) << ": calling ReadsToPathsCoreY" << endl;
               const int threads_to_use = 2 * omp_get_max_threads( ) / 3;
               ReadsToPathsCoreY( unibases, K, paths, pathsrc, 
                    pathsdb, run_dir + "/" + PATCHDIR, threads_to_use );

               cout << Date( ) << ": calling Unipath" << endl;
               Unipath( paths, pathsrc, pathsdb, unipaths, unipathsdb );
               digraph A;
               BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, 
                    unipathsdb, A );
               HyperKmerPath h;
               BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
               KmerBaseBroker kbb( K, paths, pathsrc, pathsdb, unibases );
               unibases.clear( );
               for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
                    unibases.push_back( kbb.Seq( h.EdgeObject(i) ) );
               cout << Date( ) << ": " << unibases.size( ) << " new unibases" 
                    << endl;
               cout << Date( ) << ": " << TimeSince(rebuild_clock)
                    << " rebuilding unipaths" << endl;
     
               // Write new files.

               cout << Date( ) << ": writing new files" << endl;
               unibases.WriteAll( run_dir + "/" + OUT_HEAD + ".unibases.k" + KS );
               paths.WriteAll( run_dir + "/" + OUT_HEAD + ".paths.k" + KS );
               pathsrc.WriteAll( run_dir + "/" + OUT_HEAD + ".paths_rc.k" + KS );
               unipaths.WriteAll( run_dir + "/" + OUT_HEAD + ".unipaths.k" + KS );
               BinaryWrite3( run_dir + "/" + OUT_HEAD + ".pathsdb.k" + KS, pathsdb );
               BinaryWrite3( run_dir + "/" + OUT_HEAD + ".unipathsdb.k" + KS, 
                    unipathsdb );    }    }

     // Finish up.

     cout << "\nSUMMARY:\n";
     cout << attempted_joins << " joins attempted" << endl;
     cout << Sum(joined) << " joins made, yielding " << npatches
          << " patches" << endl;
     if (log_aligns) 
          cout << Sum(joined) - Sum(perfect) << " imperfect joins made" << endl;
     cout << "\n" << Date( ) << ": done, time used = " 
          << TimeSince(clock) << endl;    }
