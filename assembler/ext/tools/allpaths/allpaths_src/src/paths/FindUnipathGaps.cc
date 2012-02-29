///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FindUnipathGaps.  Code based on RemodelGaps to estimate gaps between unipaths.
//
// There are some problems with the data structures, arising from the "placements"
// data structure, which makes it impossible to simultaneously have an alignment
// of a read to both a unipath and its reverse complement.  We've compensated for
// this in the code here in several deleterious ways:
// (1) We don't consider links from a unipath to itself or its reverse complement.
// (2) Consequently we only link between copy number one unipaths.
// (3) We copy genome-wide data structures for each gap, which is hugely inefficient.
//
// This code does not use long links, but should.
// We could try using fragment pair links too, although this might not help.
//
// This only works on haploids at present. 
//
// Note that in general, estimation of gaps should use not just adjacent unipaths,
// but also those that are farther away.
//
// Currently this code reads in preexisting gap values and compares to them.
// This is presumably temporary.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include <unistd.h>
#include "Basevector.h"
#include "FastIfstream.h"
#include "FeudalMimic.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "ParseSet.h"
#include "Set.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/PdfEntry.h"
#include "paths/RemodelGapTools.h"
#include "paths/Ulink.h"
#include "paths/UnibaseUtils.h"
#include "system/HostName.h"


// Switch the alignment of the reads from m1 to m2 ( the rc of m1 )
void FlipUniAlign( int m1, int m2,
          const vecbasevector& tigs,
          const PairsManager& pairs,
          const vec< unsigned short int >& base_lens,
          const vec< vec<longlong> > & aligns_index,
          // output
          vec<Bool> & placed_fw, 
          vec<Bool> & placed_rc,
          vec< pair<int,int> > & placement,
          // switch
          bool FLIP_HALF=false     // only flip the first read of the pair, to allow
                                   // link from unibases to the inversion of itself
     )
{
     const vec<longlong>& pids_on_tig = aligns_index[m1];
     for ( size_t i = 0; i < pids_on_tig.size(); i++ ) {
          longlong pid = pids_on_tig[i];
          if ( i >= 1 ) ForceAssertNe( pid, pids_on_tig[i-1] ); // make sure no duplication in pair ids to avoid flipping same read twice
          longlong rid1 = pairs.ID1(pid);
          longlong rid2 = pairs.ID2(pid);
          if ( ( placed_fw[rid1] || placed_rc[rid1] ) && placement[rid1].first == m1 ) {
               swap( placed_fw[rid1], placed_rc[rid1] );
               placement[rid1].second = tigs[m1].size() - base_lens[rid1] - placement[rid1].second;
               placement[rid1].first = m2;
          }
          if ( FLIP_HALF ) continue;
          if ( ( placed_fw[rid2] || placed_rc[rid2] ) && placement[rid2].first  == m1 ) {
               swap( placed_fw[rid2], placed_rc[rid2] );
               placement[rid2].second = tigs[m1].size() - base_lens[rid2] - placement[rid2].second;
               placement[rid2].first = m2;
          }
     }
}

template<int K> void FindUnipathGaps(
     const vecbasevector& unibases,
     const vecbasevector& unibases0,
     const vec<int>& to_rc,
     vec<superb>& scaffolds,
     const vec< pair<int,int> >& gaps,
     const int VERBOSITY,
     const Bool TIME_STAMPS,
     const Bool VALIDATE,
     const Bool SHOW_START_STOP,
     const String& TIGS,
     const vec<int>& tigs_to_process,
     const vecbasevector& bases,
     const PairsManager& pairs,
     const vec<unsigned short>& read_lengths,
     vec<Bool>& placed_fw,
     vec<Bool>& placed_rc,
     vec< pair<int,int> >& placement,
     const vec< vec<longlong> >& aligns_index,
     const vec<int>& true_gap,
     const vec<Bool>& true_gap_computed,
     const vec<int>& max_dist,
     const int nlibs_jump,
     const vec<int>& libs_to_use,
     const vec<Bool>& lfail,
     const vec< vec<int> >& DISTS,
     const vec< vec<double> >& X,
     const String& run_dir,
     const String& HEAD,
     vec<int>& new_gap, 
     vec<int>& new_dev,
     vec<Bool>& new_gap_computed,
     const Bool WRITE,
     const String& KS1,
     const vec< vec< pair<int,int> > >& links_to )
{
     // Build a map of kmers in the unibases.

     if (TIME_STAMPS) cout << Date( ) << ": building kmer maps for unibases" << endl;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup( unibases, kmers_plus );
     if (TIME_STAMPS) cout << Date( ) << ": copying" << endl;
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Go through the scaffolds.

     if (TIME_STAMPS) cout << Date( ) << ": going through scaffolds" << endl;
     if ( VERBOSITY >= 1 ) cout << "\n";

     // separate the gaps into different groups so that the gaps within the same
     // group don't share any unibases, which allows parallelization
     // In pass = 0, we deal with gaps from unibase to the inversion of itself,
     //              don't forget to to flip back the alignment after processing
     // In pass = 1, we deal with gaps without flipping alignment
     // In pass >= 2, we process gaps involving alignment flipping, flipped alignment
     //               will be kept to save processing time.
     //
     vec<int> gaps_to_process( gaps.size(), vec<int>::IDENTITY );
     vec<bool> aligned( unibases.size(), False );
     for ( int m = 0; m < aligned.isize(); m++ )
          if ( unibases0[m].size() > 0 ) aligned[m] = true;
     int pass = -1;
     while ( gaps_to_process.size() > 0 ) {
          pass++;
          vec<int>  gaps_in_group;
          set<int>  tigs_in_group; // list of conflicting unibases
          vec< Bool > to_del( gaps_to_process.size(), False);
          for ( int i = 0; i < gaps_to_process.isize(); i++ ) {
               int gi = gaps_to_process[i];
               int s = gaps[gi].first, jg = gaps[gi].second;
               superb& S = scaffolds[s];
               int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
               if ( TIGS != "" && !BinMember( tigs_to_process, m1 ) ) continue;
               Bool using_m1 = unibases0[m1].size( ) > 0;
               Bool using_m2 = unibases0[m2].size( ) > 0;
               if ( pass == 0 ) {
                    // In pass 1, we handle gaps from u to u*
                    if ( m1 == to_rc[m2] ) {
                         gaps_in_group.push_back( gi );
                         to_del[i] = true;
                    }
               } else if ( pass == 1 ) {
                    // In pass 1, we handle gaps that does not require alignment flipping
                    if ( using_m1 && using_m2 ) {
                         gaps_in_group.push_back( gi );
                         to_del[i] = true;
                    }
               } else {
                    if ( tigs_in_group.find(m1) != tigs_in_group.end() ) continue;
                    if ( tigs_in_group.find(m2) != tigs_in_group.end() ) continue;
                    tigs_in_group.insert( to_rc[m1] );
                    tigs_in_group.insert( to_rc[m2] );
                    // if flipping required for m1 or m2, do not allow other gaps in the same group
                    // that involves m1 or m2
                    if ( ! aligned[m1] ) tigs_in_group.insert( m1 );
                    if ( ! aligned[m2] ) tigs_in_group.insert( m2 );
                    gaps_in_group.push_back( gi );
                    to_del[i] = true;
               }
          }
          EraseIf(gaps_to_process, to_del);
          if ( VERBOSITY >= 2) {
               cout << "pass " << pass << ", " << gaps_in_group.isize() << " gaps: ";
               for ( int j = 0; j < gaps_in_group.isize(); j++ ) {
                    int gi = gaps_in_group[j];
                    int s = gaps[gi].first, jg = gaps[gi].second;
                    superb& S = scaffolds[s];
                    int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
                    cout <<  m1 << "_" << m2 << "(" << to_rc[m2] << "_" << to_rc[m1] <<  ")";
               }
               cout << endl;
               cout << "gaps left : " << gaps_to_process.size() << endl;
          }

          if ( pass == 0 ) 
               for ( int i = 0; i < gaps_in_group.isize( ); i++ ) {
                    int gi = gaps_in_group[i];
                    int s = gaps[gi].first, jg = gaps[gi].second;
                    superb& S = scaffolds[s];
                    int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
                    ForceAssert( m1 == to_rc[m2] );
                    // flip half of the aligns from m1 to m1*
                    if ( unibases0[m1].size( ) > 0 )
                         FlipUniAlign( m1, to_rc[m1], unibases, pairs, read_lengths, aligns_index, 
                                   placed_fw, placed_rc, placement, true);
               }

     #pragma omp parallel for schedule( dynamic, 1 )
     for ( int i = 0; i < gaps_in_group.isize( ); i++ )
     {    int gi = gaps_in_group[i];
          int s = gaps[gi].first, jg = gaps[gi].second;
          superb& S = scaffolds[s];
          int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
          if ( TIGS != "" && !BinMember( tigs_to_process, m1 ) ) continue;
          int tgap = 0;
          if ( VALIDATE && true_gap_computed[gi] ) tgap = true_gap[gi];
          if (SHOW_START_STOP)
          {    
               #pragma omp critical
               {    cout << "start " << m1 << endl;    }    }

          // Compute overlaps.

          vec<int> accepted_overlaps;
          int max_overlap;
          ComputeOverlaps( unibases[m1], unibases[m2], m1, max_dist,
               kmers_plus, kmers, VERBOSITY, accepted_overlaps, max_overlap );

          Bool using_m1 = unibases0[m1].size( ) > 0;
          Bool using_m2 = unibases0[m2].size( ) > 0;
          if ( pass == 1 ) ForceAssert( using_m2 && using_m1 );

          // Flip, horrible.  Needed because we only aligned to half the unibaess.
          if ( pass > 1 ) {
               if ( ! aligned[m1] ) {
                    FlipUniAlign( to_rc[m1], m1, unibases, pairs, read_lengths, aligns_index, placed_fw, placed_rc, placement);
                    aligned[m1] = true;
                    aligned[ to_rc[m1] ] = false;
               }
               if ( ! aligned[m2] ) {
                    FlipUniAlign( to_rc[m2], m2, unibases, pairs, read_lengths, aligns_index, placed_fw, placed_rc, placement );
                    aligned[m2] = true;
                    aligned[ to_rc[m2] ] = false;
               } 
          }
          // Compute gap.

          int gap;
          double dev;
          Bool gap_predicted;
          PredictGap( "jump", kmers_plus, kmers, max_overlap, nlibs_jump, 0, 
               libs_to_use, lfail, unibases, DISTS, X, bases, pairs, read_lengths, 
               aligns_index, placed_fw, placed_rc, placement, 
               gap_id( gap_id::BETWEEN, m1, m2 ), 
               VERBOSITY, gap_predicted, gap, dev );


          if ( accepted_overlaps.solo( ) && gap_predicted )
          {    int ogap = -accepted_overlaps[0];
               if ( Abs( gap - ogap ) <= 3.0 * dev ) gap = ogap;    }
          #pragma omp critical
          {    if ( VERBOSITY >= 2 ) cout << "\n";
               double offby = double( Abs( tgap - gap ) ) / dev;
               //cout << "gi " << gi << endl;
               if ( VERBOSITY >= 1 )
               {    cout << "m1 = " << m1 << ", m2 = " << m2;
                    if (gap_predicted)
                    {    cout << ", gap = " << gap << " +/- "
                              << setiosflags(ios::fixed) << setprecision(1) 
                              << dev;    }
                    else cout << ", failed to predict gap";
                    if ( VALIDATE && true_gap_computed[gi] )
                    {    cout << ", tgap = " << tgap;
                         if (gap_predicted) cout << ", offby = " << offby;    }
                    if ( S.Gap(jg) == 999 && S.Dev(jg) == 99 )
                         cout << ", OLD MISSING";
                    cout << endl;    }    }
          if (gap_predicted)
          {    new_gap[gi] = gap, new_dev[gi] = int(round(dev));    }
          new_gap_computed[gi] = gap_predicted;    } 

          // flip back the u to u* aligns
          if ( pass == 0 ) 
               for ( int i = 0; i < gaps_in_group.isize( ); i++ ) {
                    int gi = gaps_in_group[i];
                    int s = gaps[gi].first, jg = gaps[gi].second;
                    superb& S = scaffolds[s];
                    int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
                    ForceAssert( m1 == to_rc[m2] );
                    // flip back
                    if ( unibases0[m1].size( ) > 0 )
                         FlipUniAlign( to_rc[m1], m1, unibases, pairs, read_lengths, aligns_index, 
                                   placed_fw, placed_rc, placement, true);
               }
     }
     if (TIME_STAMPS) cout << Date( ) << ": Finished gap calculation." << endl;

     // Generate report and write results.

     vec<int> devs, adevs;
     vec<double> offbys, aoffbys;
     vec< triple<int,int,int> > results;
     int fails = 0;
     for ( int gi = 0; gi < gaps.isize( ); gi++ )
     {    int s = gaps[gi].first, jg = gaps[gi].second;
          superb& S = scaffolds[s];
          int m1 = S.Tig(jg), m2 = S.Tig(jg+1);
          if ( TIGS != "" && !BinMember( tigs_to_process, m1 ) ) continue;
          if ( VALIDATE && true_gap_computed[gi] )
          {    if ( !new_gap_computed[gi] )
               {    fails++;
                    continue;    }
               int gap = new_gap[gi], dev = new_dev[gi], tgap = true_gap[gi];
               results.push( gap, tgap, S.Gap(jg) );
               devs.push_back(dev);
               adevs.push_back( S.Dev(jg) );
               double aoffby = double( Abs( S.Gap(jg) - tgap ) ) / double(S.Dev(jg));
               double offby = double( Abs( tgap - gap ) ) / dev;
               offbys.push_back(offby), aoffbys.push_back(aoffby);    }    }
     for ( int gi = 0; gi < gaps.isize( ); gi++ )
     {    if ( new_gap_computed[gi] )
          {    int s = gaps[gi].first, jg = gaps[gi].second;
               superb& S = scaffolds[s];
               S.SetGap( jg, new_gap[gi] ), S.SetDev( jg, new_dev[gi] );    }    }
     if ( VALIDATE && ( TIGS == "" || tigs_to_process.size( ) > 1 ) )
     {    WriteReportAboutGaps( devs, adevs, offbys, aoffbys, results, true_gap,
               true_gap_computed, fails, tigs_to_process );    }

     int nself = 0, nrc = 0, ntotal=0;
     if (WRITE)
     {    Ofstream( out, run_dir + "/" + HEAD + ".unibases.k" + KS1 
               + ".predicted_gaps.txt" );
          out << "#u1 u2 sep dev nlinks " << endl;
          for ( int gi = 0; gi < gaps.isize( ); gi++ )
          {    int s = gaps[gi].first, jg = gaps[gi].second;
               superb& S = scaffolds[s];
               int u1 = S.Tig(jg), u2 = S.Tig(jg+1);
               if ( new_gap_computed[gi] )
               {    int sep = new_gap[gi], dev = new_dev[gi];
                    ntotal++;
                    if ( u1 == u2 ) { nself++; cout << u1 <<"_"<<u1 << endl; }
                    if ( u1 == to_rc[u2] ) nrc++;
                    int nlinks = -1;
                    for ( int j = 0; j < links_to[u1].isize( ); j++ )
                    {    if ( links_to[u1][j].first == u2 )
                         {    nlinks = links_to[u1][j].second;
                              break;    }    }
                    out << u1 << " " << u2 << " " << sep << " " << dev
                         << " " << nlinks << "\n";    }    }    }
             
     cout <<  endl;
     cout <<  ntotal << " gaps are calculated. " << endl;
     cout <<  nself << " gaps are self-linking " << endl;
     cout <<  nrc << " gaps are between unibases and their self-inversions" 
          << endl;    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(HEAD, "extended40.shaved.patched");
     CommandArgument_Int_OrDefault(K1, 96);
     CommandArgument_Int_OrDefault(K, 40);
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault_Doc(TIGS, "",
          "if specified, only process gaps have these contigs on their left");
     CommandArgument_String_OrDefault(FRAG_READS_IN, "frag_reads_filt");
     CommandArgument_String_OrDefault(JUMP_READS_IN, "jump_reads_ec");
     CommandArgument_String_OrDefault(LONG_JUMP_READS_IN, "long_jump_reads_filt");
     CommandArgument_String_OrDefault_Doc(LIB_TYPES, "jump",
          "theoretically could be {frag,jump,long}, see comments in code");
     CommandArgument_String_OrDefault(LIBS, "");
     CommandArgument_Bool_OrDefault(WRITE, True);
      CommandArgument_Int_OrDefault_Doc( MIN_KMERS, 800,
	"Min number of k-mers in the unipath used in scaffolding." );
     // Logging

     CommandArgument_Int_OrDefault_Doc(VERBOSITY, -1,
          "VERBOSITY=1: print one line for each gap\n"
          "VERBOSITY=2: print more data for each gap\n"
          "VERBOSITY=3: print a lot more\n"
          "VERBOSITY=4: print even more\n"
          "(Note VERBOSITY>=2 should only be used if TIGS is set to one contig.)");
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(SHOW_START_STOP, False);
     CommandArgument_Bool_OrDefault(TIME_STAMPS, True);
     CommandArgument_Bool_OrDefault(HISTOGRAMS, False);
     CommandArgument_Bool_OrDefault(PRINT_LINKS, False);

     EndCommandArguments;

     // Define directories, etc.

     if (TIME_STAMPS) cout << Date( ) << ": begin" << endl;
     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     vec<String> lib_types_to_use;
     ParseStringSet( LIB_TYPES, lib_types_to_use );
     Bool USE_FRAGS = Member( lib_types_to_use, String("frag") );
     vec<int> libs_to_use;
     ParseIntSet( LIBS, libs_to_use );

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Heuristics.

     const int min_kmers = MIN_KMERS;
     const int min_kmers1_self = 25000;
     const int min_links = 2;
     const int min_links_self = 10;

     // Define unibases to process.

     vec<int> tigs_to_process;
     if ( TIGS != "" ) 
     {    ParseIntSet( TIGS, tigs_to_process );
          if ( VERBOSITY < 0 ) VERBOSITY = 1;    }
     if ( VERBOSITY >= 2 && !tigs_to_process.solo( ) )
     {    cout << "If you use VERBOSITY >= 2, then you need to specify "
               << "a single contig.\n";
          exit(1);    }

     // Load unibases and ancillary data structures.

     String KS1 = ToString(K1);
     vecbasevector unibases( run_dir + "/" + HEAD + ".unibases.k" + KS1 );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc, K1 );

     // Load predicted gaps between unibases.

     vec<ulink_with_uids> ulinks;
     vec< vec<int> > ulinks_index;
     if (VALIDATE)
     {    ulinks_index.resize( unibases.size( ) );
          String line;
          fast_ifstream in( 
               run_dir + "/" + HEAD + ".unipaths.predicted_gap.k" + KS1 );
          while(1)
          {    getline( in, line );
               if ( in.fail( ) ) break;
               if ( line.Contains( "#", 0 ) ) continue;
               int u1, u2, sep, dev, nlinks;
               istringstream iline( line.c_str( ) );
               iline >> u1 >> u2 >> sep >> dev >> nlinks;
               int nkmers1 = unibases[u1].isize( ) - K1 + 1;
               int nkmers2 = unibases[u2].isize( ) - K1 + 1;
               ulinks_index[u1].push_back( ulinks.size( ) );
               ulinks.push( u1, u2, sep, dev, 0, 0, nlinks );    }    }
     
     // Align the jumping reads to these unibases.  We first create a copy of the
     // unibases in which reverse complements are deleted.  Otherwise, since we
     // require unique placement, nothing would align.

     vecbasevector unibases0(unibases);
     for ( int u = 0; u < (int) unibases.size( ); u++ )
     {    if ( to_rc[u] < u ) unibases0[u].resize(0);
               }
     String jhead = JUMP_READS_IN;
     PairsManager pairs( run_dir + "/" + jhead + ".pairs" );
     vecbasevector bases( run_dir + "/" + jhead + ".fastb" );
     vec<Bool> placed_fw, placed_rc; 
     vec< pair<int,int> > placement;
     vec< vec<longlong> > aligns_index;
     AlignReads( K, 20, True, unibases0, bases, pairs, placed_fw, placed_rc, 
          placement, aligns_index, False );
     for ( int u = 0; u < (int) unibases.size( ); u++ )
          if ( to_rc[u] < u ) aligns_index[u] = aligns_index[ to_rc[u] ];

     // Count number of libraries.

     int nlibs_frag, nlibs_jump, nlibs_long;
     GetLibCounts( run_dir, lib_types_to_use, 0, False, FRAG_READS_IN, 
          JUMP_READS_IN, LONG_JUMP_READS_IN, nlibs_frag, nlibs_jump, nlibs_long );
     int nlibs = nlibs_frag + nlibs_jump + nlibs_long;

     // Define distributions for the libraries.

     vec<unsigned short> read_lengths( bases.size( ) );
     for ( size_t i = 0; i < bases.size( ); i++ )
          read_lengths[i] = bases[i].size( );
     vec<Bool> lfail;
     vec<int> max_dist;
     vec< vec<int> > DISTS;
     vec< vec<double> > X;
     DefineDistributions( nlibs_jump, 0, unibases0, read_lengths, pairs, placed_fw, 
          placed_rc, placement, aligns_index, TIME_STAMPS, HISTOGRAMS, lfail, 
          max_dist, DISTS, X );

     // Define linked unipaths.

     vec< vec< pair<int,int> > > links_to( unibases.size( ) );
     for ( int u1 = 0; u1 < (int) unibases.size( ); u1++ )
     {    if ( to_rc[u1] < u1 ) continue;
          int nkmers1 = unibases[u1].isize( ) - K1 + 1;
          if ( nkmers1 < min_kmers ) continue;
          vec< pair<int,int> > linked_fw, linked_rc;
          // !!! Note that aligns_index contains all pairs that are aligned to u1.
          // If boths reads are aligned to u1, the pair id will appear twice here,
          // and can possibly cause double counting if not handled properly.
          for ( int j = 0; j < aligns_index[u1].isize( ); j++ )
          {    longlong pid = aligns_index[u1][j];
               longlong id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);

	       // force id1 to be the read aligned to u1
	       if ( placement[id1].first != u1 ) swap( id1, id2 );
	       ForceAssertEq( u1, placement[id1].first );

               int pos1 = placement[id1].second;
               if ( !placed_fw[id2] && !placed_rc[id2] ) continue;
               int u2 = placement[id2].first, pos2 = placement[id2].second;

               // We don't allow links from a unibase to itself or its reverse
               // complement.  Note that this is OK for CN1 unibases but not in
               // general.
               //  !! Now we allow both u to u and u to u* links

               // if ( /* u2 == u1 || */ u2 == to_rc[u1] ) continue;
               // u2 == to_rc[u1] shall not happen as only one copy of each unibases are used.
               ForceAssertNe( u2, to_rc[u1] );

               if ( u2 == u1 )
               {    //if ( placement[id1].first != u1 ) continue;
                    if ( nkmers1 < min_kmers1_self ) continue;    }

               Bool fw1 = placed_fw[id1], fw2 = placed_fw[id2];
               int nkmers2 = unibases[u2].isize( ) - K1 + 1;
               if ( nkmers2 < min_kmers ) continue;
               int libid = pairs.libraryID(pid);
               int dist1 = ( fw1 ? unibases[u1].isize( ) - pos1 : pos1 );
               int dist2 = ( fw2 ? unibases[u2].isize( ) - pos2 : pos2 );
               if ( dist1 > max_dist[libid] ) continue;
               if ( dist2 > max_dist[libid] ) continue;
               int antidist1 = ( !fw1 ? unibases[u1].isize( ) - pos1 : pos1 );
               int antidist2 = ( !fw2 ? unibases[u2].isize( ) - pos2 : pos2 );
               int antidist = antidist1 + antidist2;
               if (PRINT_LINKS) 
               {    int v1 = ( fw1 ? u1 : to_rc[u1] );
                    int v2 = ( !fw2 ? u2 : to_rc[u2] );
                    cout << "LINK: ";
                    PRINT8( v1, v2, dist1, dist2, antidist, id1, id2, pid );    }
               ( fw1 ? linked_fw : linked_rc )
                    .push( !fw2 ? u2 : to_rc[u2], antidist );    }
          Sort(linked_fw), Sort(linked_rc);
          // Parameter to avoid innie links:
          const int min_antidist = 500; // should not be hardcoded!!!!!!!!!!!!!!!!!!
          for ( int i = 0; i < linked_fw.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < linked_fw.isize( ); j++ )
                    if ( linked_fw[j].first != linked_fw[i].first ) break;
               if ( j - i >= min_links 
                    && ( u1 != linked_fw[i].first || j - i >= min_links_self ) ) 
               {    vec<int> antidists;
                    for ( int k = i; k < j; k++ )
                         antidists.push_back( linked_fw[k].second );
                    Sort(antidists);
                    if ( Median(antidists) >= min_antidist )
                         links_to[u1].push( linked_fw[i].first, j-i );    }
               i = j - 1;    }
          for ( int i = 0; i < linked_rc.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < linked_rc.isize( ); j++ )
                    if ( linked_rc[j].first != linked_rc[i].first ) break;
               if ( j - i >= min_links && ( to_rc[u1] != linked_rc[i].first 
                    || j - i >= min_links_self ) )
               {    vec<int> antidists;
                    for ( int k = i; k < j; k++ )
                         antidists.push_back( linked_rc[k].second );
                    Sort(antidists);
                    if ( Median(antidists) >= min_antidist )
                         links_to[ to_rc[u1] ].push( linked_rc[i].first, j-i );    }
               i = j - 1;    }    }
     for ( int u1 = 0; u1 < (int) unibases.size( ); u1++ )
     {    vec< pair<int,int> >& L = links_to[u1], L2;
          Sort(L);
          for ( int j = 0; j < L.isize( ); j++ )
          {    int k;
               for ( k = j + 1; k < L.isize( ); k++ )
                    if ( L[k].first != L[j].first ) break;
               int nlinks = 0;
               for ( int r = j; r < k; r++ )
                    nlinks += L[r].second;
               if ( nlinks >= min_links ) L2.push( L[j].first, nlinks );
               j = k - 1;    }
          L = L2;    }

     // Make goofy scaffolds representing these links.

     vec<superb> scaffolds;
     for ( int u1 = 0; u1 < (int) unibases.size( ); u1++ )
     {    for ( int i = 0; i < links_to[u1].isize( ); i++ )
          {    int u2 = links_to[u1][i].first;
               superb s;
               s.SetNtigs(2);
               s.SetTig( 0, u1 );
               s.SetTig( 1, u2 );
               int sep = 999, dev = 99; // default, means no value
               if (VALIDATE)
               {    for ( int j = 0; j < ulinks_index[u1].isize( ); j++ )
                    {    const ulink_with_uids& l = ulinks[ ulinks_index[u1][j] ];
                         if ( l.u2 == u2 )
                         {    sep = l.sep, dev = l.dev;    }    }    }
               s.SetGap(0, sep), s.SetDev(0, dev);
               scaffolds.push_back(s);    }    }

     // Load genome and find true gaps.  We only use the 1000-mers flanking the 
     // gaps have unique placements in the reference.  We seed on 100-mers to 
     // speed up the computation.

     vecbasevector genome;
     vec<int> true_gap;
     vec<Bool> true_gap_computed;
     if (VALIDATE) 
     {    genome.ReadAll( data_dir + "/genome.fastb" );
          FindTrueGaps( genome, unibases, scaffolds, tigs_to_process, TIME_STAMPS,
               true_gap, true_gap_computed );    }

     // Define gaps.

     vec< pair<int,int> > gaps;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int jg = 0; jg < S.Ngaps( ); jg++ )
               gaps.push( s, jg );    }
     vec<int> new_gap( gaps.size( ) ), new_dev( gaps.size( ) );
     vec<Bool> new_gap_computed( gaps.size( ), False );

     // Run the rest of the program.

     if ( K == 40 )
     {    FindUnipathGaps<40>( unibases, unibases0, to_rc, scaffolds, gaps,
               VERBOSITY, TIME_STAMPS, VALIDATE, SHOW_START_STOP, TIGS,
               tigs_to_process, bases, pairs, read_lengths, placed_fw, placed_rc,
               placement, aligns_index, true_gap, true_gap_computed, max_dist,
               nlibs_jump, libs_to_use, lfail, DISTS, X, run_dir, HEAD,
               new_gap, new_dev, new_gap_computed, WRITE, KS1, links_to );    }
     else if ( K == 80 )
     {    FindUnipathGaps<80>( unibases, unibases0, to_rc, scaffolds, gaps,
               VERBOSITY, TIME_STAMPS, VALIDATE, SHOW_START_STOP, TIGS,
               tigs_to_process, bases, pairs, read_lengths, placed_fw, placed_rc,
               placement, aligns_index, true_gap, true_gap_computed, max_dist,
               nlibs_jump, libs_to_use, lfail, DISTS, X, run_dir, HEAD,
               new_gap, new_dev, new_gap_computed, WRITE, KS1, links_to );    }
     else 
     {    cout << "FindUnipathGaps not implemented for K = " << K << endl;
          exit(1);    }
     cout << "\ntime used = " << TimeSince(clock) << " (" << getHostName()
             << ")\n\n";
     cout << "command: " << command.TheCommand( ) << "\n\n";    }
