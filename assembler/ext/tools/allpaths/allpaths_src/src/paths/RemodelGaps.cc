///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// RemodelGaps.  Estimate gaps between contigs in scaffolds.
//
// Known problems:
//
// 1. Some of the algorithm is hardcoded to use only jumps.
// 2. Not tested on large genomes.

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
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/RemodelGapTools.h"
#include "system/HostName.h"

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".remodel");
     CommandArgument_Int_OrDefault_Doc(NUM_THREADS, -1,
          "number of threads to be used (use all available if negative)");
     CommandArgument_String_OrDefault_Doc(TIGS, "",
          "if specified, only process gaps have these contigs on their left");
     CommandArgument_String_OrDefault(FRAG_READS_IN, "frag_reads_filt");
     CommandArgument_String_OrDefault(JUMP_READS_IN, "jump_reads_filt");
     CommandArgument_String_OrDefault(LONG_JUMP_READS_IN, "long_jump_reads_filt");
     CommandArgument_String_OrDefault_Doc(LIB_TYPES, "jump",
          "theoretically could be {frag,jump,long}, see comments in code");
     CommandArgument_String_OrDefault(LIBS, "");
     CommandArgument_Bool_OrDefault(WRITE, False);

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

     EndCommandArguments;

     // Define directories, etc.

     if (TIME_STAMPS) cout << Date( ) << ": begin" << endl;
     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String head = sub_dir + "/" + SCAFFOLDS_IN;
     vec<String> lib_types_to_use;
     ParseStringSet( LIB_TYPES, lib_types_to_use );
     Bool USE_FRAGS = Member( lib_types_to_use, String("frag") );
     vec<int> libs_to_use;
     ParseIntSet( LIBS, libs_to_use );

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load assembly.

     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     vecbasevector tigs( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     int ntigs = tigs.size( );
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     vec<int> tigs_to_process;
     if ( TIGS != "" ) 
     {    ParseIntSet( TIGS, tigs_to_process );
          if ( VERBOSITY < 0 ) VERBOSITY = 1;    }
     if ( VERBOSITY >= 2 && !tigs_to_process.solo( ) )
     {    cout << "If you use VERBOSITY >= 2, then you need to specify "
               << "a single contig.\n";
          exit(1);    }

     // Load genome and find true gaps.  We only use the 1000-mers flanking the 
     // gaps have unique placements in the reference.  We seed on 100-mers to 
     // speed up the computation.

     vecbasevector genome;
     vec<int> true_gap;
     vec<Bool> true_gap_computed;
     if (VALIDATE) 
     {    genome.ReadAll( data_dir + "/genome.fastb" );
          FindTrueGaps( genome, tigs, scaffolds, tigs_to_process, TIME_STAMPS,
               true_gap, true_gap_computed );    }

     // Build a map of kmers in the contigs.

     if (TIME_STAMPS) cout << Date( ) << ": building kmer maps for contigs" << endl;
     const int K = 40;
     vec< triple<kmer<K>,int,int> > kmers_plus;
     MakeKmerLookup( tigs, kmers_plus );
     if (TIME_STAMPS) cout << Date( ) << ": copying" << endl;
     vec< kmer<K> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Count number of libraries.

     int nlibs_frag, nlibs_jump, nlibs_long;
     GetLibCounts( run_dir, lib_types_to_use, VERBOSITY, TIME_STAMPS,
		   FRAG_READS_IN, JUMP_READS_IN, LONG_JUMP_READS_IN,
		   nlibs_frag, nlibs_jump, nlibs_long );
     int nlibs = nlibs_frag + nlibs_jump + nlibs_long;

     // Define gaps.

     vec< pair<int,int> > gaps;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int jg = 0; jg < S.Ngaps( ); jg++ )
               gaps.push( s, jg );    }
     vec<int> new_gap( gaps.size( ) ), new_dev( gaps.size( ) );
     vec<Bool> new_gap_computed( gaps.size( ), False );

     // Go through two passes, one for fragment reads and one for jumps.

     for ( int pass = 0; pass < lib_types_to_use.isize( ); pass++ )
     {    if (TIME_STAMPS) cout << Date( ) << ": starting pass " << pass+1 << endl;

          // Align reads.

          if (TIME_STAMPS) cout << Date( ) << ": loading reads" << endl;
          String head;
          if ( lib_types_to_use[pass] == "frag" ) head = FRAG_READS_IN;
          if ( lib_types_to_use[pass] == "jump" ) head = JUMP_READS_IN;
          PairsManager pairs( run_dir + "/" + head + ".pairs" );
          vecbasevector bases;
          if ( lib_types_to_use[pass] == "frag" )
               bases.ReadAll( run_dir + "/" + head + ".fastb" );
          vec<unsigned short> read_lengths;
          vec<Bool> placed_fw, placed_rc; 
          vec< pair<int,int> > placement;
          vec< vec<longlong> > aligns_index;
          {    vecbasevector bases( run_dir + "/" + head + ".fastb" );
               if ( lib_types_to_use[pass] == "frag" )
               {    vec< triple<kmer<K>,int,int> > bases_kmers_plus;
                    MakeKmerLookup( bases, bases_kmers_plus );
                    vec< kmer<K> > bases_kmers( bases_kmers_plus.size( ) );
                    for ( size_t i = 0; i < bases_kmers.size( ); i++ )
                         bases_kmers[i] = bases_kmers_plus[i].first;
                    kmer<K> x, xrc;
                    for ( size_t i = 0; i < bases.size( ); i++ )
                    {    if ( bases[i].isize( ) >= K )
                         {    x.SetToSubOf( bases[i], 0 );
                              xrc = x;
                              xrc.ReverseComplement( );
                              if ( xrc < x ) x = xrc;
                              int64_t low = LowerBound( bases_kmers, x );
                              int64_t high = UpperBound( bases_kmers, x );
                              if ( high - low == 1 ) 
                                   bases[i].resize(0);    }    }    }
               read_lengths.resize( bases.size( ) );
               for ( size_t i = 0; i < bases.size( ); i++ )
                    read_lengths[i] = bases[i].size( );
               if (TIME_STAMPS) cout << Date( ) << ": aligning reads" << endl;
               AlignReads( 40, 20, lib_types_to_use[pass] == "jump", 
                    tigs, bases, pairs, placed_fw, placed_rc, 
                    placement, aligns_index, TIME_STAMPS );    }

          // Define distributions for libraries.

          vec<Bool> lfail;
          vec<int> max_dist;
          vec< vec<int> > DISTS;
          vec< vec<double> > X;
          // int lib_offset = ( pass == 2 ? nlibs_frag : 0 );
          DefineDistributions( 
               lib_types_to_use[pass] == "frag" ? nlibs_frag : nlibs_jump,
               0, tigs, 
               read_lengths, pairs, placed_fw, placed_rc, placement, aligns_index, 
               TIME_STAMPS, HISTOGRAMS, lfail, max_dist, DISTS, X );

          // Go through the scaffolds.

          if (TIME_STAMPS) cout << Date( ) << ": going through scaffolds" << endl;
          if ( VERBOSITY >= 1 ) cout << "\n";
          #pragma omp parallel for schedule( dynamic, 1 )
          for ( int gi = 0; gi < gaps.isize( ); gi++ )
          {    int s = gaps[gi].first, jg = gaps[gi].second;
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
               ComputeOverlaps( tigs[m1], tigs[m2], m1, max_dist,
                    kmers_plus, kmers, VERBOSITY, accepted_overlaps, max_overlap );

               // Compute gap.

               int gap;
               double dev;
               Bool gap_predicted;
               PredictGap( lib_types_to_use[pass], kmers_plus, kmers, max_overlap, 
                    lib_types_to_use[pass] == "frag" ? nlibs_frag : nlibs_jump,
                    0, libs_to_use, lfail, tigs,
                    DISTS, X, bases, pairs, read_lengths, aligns_index, placed_fw, 
                    placed_rc, placement, gap_id( gap_id::BETWEEN, m1, m2 ), 
                    VERBOSITY, gap_predicted, gap, dev );
               if ( accepted_overlaps.solo( ) && gap_predicted )
               {    int ogap = -accepted_overlaps[0];
                    if ( Abs( gap - ogap ) <= 3.0 * dev ) gap = ogap;    }
               #pragma omp critical
               {    if ( VERBOSITY >= 2 ) cout << "\n";
                    double offby = double( Abs( tgap - gap ) ) / dev;
                    if ( VERBOSITY >= 1 )
                    {    cout << "using " << 
                              lib_types_to_use[pass] << "s"
                         /* ( pass == 1 ? "frags" : "jumps" ) */
                              << ", m1 = " << m1;
                         if (gap_predicted)
                         {    cout << ", gap = " << gap << " +/- "
                                   << setiosflags(ios::fixed) << setprecision(1) 
                                   << dev;    }
                         else cout << ", failed to predict gap";
                         if ( VALIDATE && true_gap_computed[gi] )
                         {    cout << ", tgap = " << tgap;
                              if (gap_predicted) cout << ", offby = " << offby;    }
                         cout << endl;    }    }
               if (gap_predicted)
               {    new_gap[gi] = gap, new_dev[gi] = int(round(dev));    }
               new_gap_computed[gi] = gap_predicted;    }    }

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
     if (WRITE)
     {    String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
          vec<efasta> tigse;
          LoadEfastaIntoStrings( efasta_file, tigse );
          Assembly A( scaffolds, tigse );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );    }
     cout << "\ntime used = " << TimeSince(clock) << " (" << getHostName()
             << ")\n\n";
     cout << "command: " << command.TheCommand( ) << "\n\n";    }
