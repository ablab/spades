//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongReadPostPatcher.  Patch gaps in scaffolds using long reads.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "Fastavector.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "PackAlign.h"
#include "ParseSet.h"
#include "PrintAlignment.h"
#include "Superb.h"
#include "VecUtilities.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "paths/AssemblyEdit.h"
#include "paths/LongReadTools.h"
#include "paths/LongReadPatchOptimizer.h"

static inline 
String Tag(String S = "LRPP") { return Date() + " (" + S + "): "; } 



Bool cmp_second( const triple<int,int,int>& t1, const triple<int,int,int>& t2 )
{    return t1.second < t2.second;    }

Bool cmp_third( const triple<int,int,int>& t1, const triple<int,int,int>& t2 )
{    return t1.third < t2.third;    }

void PrintGenomicPlacements( ostream& out, basevector b, const int L,
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs )
{    vec<placementx> places = FindGenomicPlacements( b, L, genome, Glocs );
     for ( int i = 0; i < places.isize( ); i++ )
     {    const placementx& p = places[i];
          out << ( p.fw ? "fw" : "rc" ) << " matches genome " << p.g << "." << p.pos
               << "-" << p.pos + b.isize( ) << " perfectly\n";    }    }










// BestGoodFirst: calculate the longest left part of an alignment that has a
// substitution rate <= a given bound.

double BestGoodFirst( const align& a, const basevector& b1, 
     const basevector& b2, const double max_err )
{    int p1 = a.pos1( ), p2 = a.pos2( ), errs = 0, best_len = 0;
     for ( int j = 0; j < a.Nblocks( ); j++ ) 
     {    if ( a.Gaps(j) > 0 ) p2 += a.Gaps(j);
          if ( a.Gaps(j) < 0 ) p1 -= a.Gaps(j);
          for ( int x = 0; x < a.Lengths(j); x++ ) 
          {    if ( b1[p1] != b2[p2] ) errs++;
               ++p1; ++p2;    
               if ( double(errs)/double(p1) <= max_err ) best_len = p1;    }    }
     return best_len;    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     // Assembly.
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	   "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault(SCAFFOLDS_IN, 
          "linear_scaffolds0.clean.patched");
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".longread");
     CommandArgument_String_OrDefault(READS, "long_reads_orig");
     // Import a pacbio data set directly.
     CommandArgument_String_OrDefault(PACBIO_RUNS, "");
     CommandArgument_String_OrDefault_Doc(CONTIG_IDS, "all",
          "list of contigs in ParseIntSet format or \"all\"");
     CommandArgument_String_OrDefault_Doc(IDS, "all",
          "list of reads in ParseIntSet format or \"all\"");
     CommandArgument_Bool_OrDefault(CHECKPOINT, False);
     // Algorithmic heuristics.
     // Evaluation.
     CommandArgument_Bool_OrDefault_Doc(VALIDATE, False,
          "assess results versus reference genome");
     CommandArgument_Bool_OrDefault_Doc(VALIDATE_PATCHES, False,
          "assess patches versus reference genome");
     // Logging.
     CommandArgument_Bool_OrDefault_Doc(ANNOUNCE, False,
          "announce begin and end for each read");
     CommandArgument_Bool_OrDefault_Doc(DIRECT, False, 
          "directly log to cout, only a good idea for single reads");
     CommandArgument_String_OrDefault_Doc(READ_LOG, "none",
          "which reads to log details for, can be all, hard, or none");
     CommandArgument_Bool_OrDefault_Doc(SHOW_ALIGNS, False, 
          "display raw alignments, which are voluminous");
     CommandArgument_Bool_OrDefault_Doc(PERF_STATS, False, 
          "print stuff relating to tracking computational performance");
     CommandArgument_Bool_OrDefault_Doc(LOG_PROGRESS, True, 
          "log progress of computation");
     CommandArgument_Bool_OrDefault_Doc(PATCH_STATS, False, 
          "print patch stats");
     CommandArgument_Int_OrDefault(PATCH_VERBOSITY, 1);
     CommandArgument_Bool_OrDefault(PRINT_ALIGNS, False);
     // Saving.
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault_Doc(DUMP_LMR, False,
          "dump left, middle, and right sequences");
     // Other.
     CommandArgument_Bool_OrDefault(PATCH, True);
     EndCommandArguments;

     // Start, define directories.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;

     const String fn_patchers = sub_dir + "/" + 
       SCAFFOLDS_IN + "-" + READS + ".patchers";

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Load scaffolds.

     vec<superb> scaffolds;
     ReadSuperbs( sub_dir + "/" + SCAFFOLDS_IN + ".superb", scaffolds );
     int total_gaps = 0;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
          total_gaps += scaffolds[s].Ngaps( );
     String outputFile = sub_dir + "/" + SCAFFOLDS_OUT + ".edits";
     if ( total_gaps == 0 )
     {    cout << "There are no gaps to patch, so we won't do any work." << endl;
          if (WRITE)
          {    cout << "\n" << Tag( ) << ": saving patches" << endl;
               BinaryWriter::writeFile( outputFile, vec<assembly_edit>( ) );    }
          cout << "0 gaps patched" << endl << endl << Tag() 
               << "Done, time used = " << TimeSince(clock) << endl;    }

     // Load data.

     vecbasevector R; 
     if (PACBIO_RUNS != "") {
       vec<int> runs;
       ParseIntSet( PACBIO_RUNS, runs );
       for ( int i = 0; i < runs.isize( ); i++ )
       {    vecbasevector x;
            String pb_pre = ToString( runs[i] );
            pb_pre.resize(2);
            pb_pre = "0" + pb_pre;
            String dir = "/seq/pacbio_results/userdata/jobs/" + pb_pre + "/0"
                + ToString( runs[i] ) + "/data";
            String fn;
            if ( IsRegularFile( dir + "/filtered_subreads.fa" ) )
                 fn = dir + "/filtered_subreads.fa";
            else fn = dir + "/filtered_subreads.fasta";
	    FetchReads( x, 0, fn );
	    R.Append(x);    }
     } else {
       R.ReadAll( run_dir + "/" + READS + ".fastb");
     }

     // Use only reads specified by command line
     vec<int> ids;
     if ( IDS == "all" )
     {    for ( size_t j = 0; j < R.size( ); j++ )
               ids.push_back(j);    }
     else {
       ParseIntSet( IDS, ids );
       for ( int j = 0; j < ids.isize( ); j++ ) {
         ForceAssertLe( ids[j], int(R.size()) );
         // cout << "ids[" << j << "]= " << ids[j] << endl;
       }
     }

     vecbasevector genome;
     if (VALIDATE || VALIDATE_PATCHES) genome.ReadAll( data_dir + "/genome.fastb" );

     // Define heuristic constants for validation.

     const int LG = 12;
     const int G_initial = 100;
     const int G_initial_bandwidth_div = 10;
     const double G_initial_error_rate_max = 0.25;
     const double G_bandwidth_frac = 0.2;
     const int validation_twiddle = 10;

     // Define heuristic constants for de novo algorithm.

     const int L = 12;
     const int flank = 19;
     const int max_errs = 10;
     double max_error_rate = 0.26;
     const double bandwidth_div = 5;

     // Load the assembly.

     vecbasevector T( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     vecbasevector Trc(T);
     ReverseComplement(Trc);
     
     int ntigs = T.size( );
     vec<int> idsT;
     if ( CONTIG_IDS == "all" ) {
       for (int i = 0; i != ntigs; i++) 
         idsT.push_back(i);
     }
     else {
       ParseIntSet( CONTIG_IDS, idsT );
       cout << "looking only at contigs: ";
       for ( int j = 0; j < idsT.isize( ); j++ )
         cout << " " << idsT[j];
       cout << endl;
     }

     String ASSEMBLY = "linear_scaffolds0.clean.patched";
     vec<int> to_super( ntigs, -1 ), to_super_pos( ntigs, -1 );
     for ( int i = 0; i < scaffolds.isize( ); i++ ) 
     {    for ( int j = 0; j < scaffolds[i].Ntigs( ); j++ ) 
          {    to_super[ scaffolds[i].Tig(j) ] = i;
               to_super_pos[ scaffolds[i].Tig(j) ] = j;    }    }
     String tigsa_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fasta";
     vec<fastavector> tigsa;
     LoadFromFastaFile( tigsa_file, tigsa );

     // Check arguments.

     ForceAssert( READ_LOG == "all" || READ_LOG == "hard" || READ_LOG == "none" );

     // Report reads stats.

     longlong total_length = 0;
     for ( size_t i = 0; i < R.size( ); i++ )
          total_length += R[i].size( );
     cout << "have " << ToStringAddCommas( R.size( ) ) 
          << " reads of mean length " << total_length / R.size( ) << "\n";
     if (VALIDATE)
     {    longlong gbases = 0;
          for ( size_t i = 0; i < genome.size( ); i++ )
               gbases += genome[i].size( );
          cout << ", and providing " << setprecision(3) 
               << double(total_length)/double(gbases) << "x coverage";    }
     cout << endl;

     // Add the reverse complements of reads.

     vecbasevector R2;
     for ( size_t i = 0; i < R.size( ); i++ )
     {    R2.push_back( R[i] );
         R[i].ReverseComplement( );
         R2.push_back( R[i] );    }
     R = R2;

     // Hash the genome.

     vec< vec< pair<int,int> > > Glocs;
     if (VALIDATE || VALIDATE_PATCHES)
     {    Glocs.resize( IPow( 4, LG ) );
          for ( size_t i = 0; i < genome.size( ); i++ )
          {    for ( int j = 0; j <= genome[i].isize( ) - LG; j++ )
               {    int n = KmerId( genome[i], LG, j );
                    Glocs[n].push( i, j );    }    }    }

     // Hash the contigs.

     vec< vec< pair<int,int> > > Tlocs( IPow( 4, L ) );
     for ( size_t i = 0; i < idsT.size( ); i++ ) {
       size_t iT = idsT[i];
       for ( int j = 0; j <= T[iT].isize( ) - L; j++ ) {
         int n = KmerId( T[iT], L, j );
         Tlocs[n].push( iT, j );
       }
     }

     // Set up global stats.

     int simple_good = 0, found_alls = 0;

     // Go through the reads.

     if (!CHECKPOINT || !IsRegularFile(fn_patchers)) {

       vec<GapPatcher> patchers;
       vec<int> best_len;
       vec<String> report( R.size( ) );

       uint reads_processed = 0;
       cout << Tag() << ": going through the reads (100 dots to follow)\n";
       #pragma omp parallel for schedule(dynamic, 1)
       for ( size_t zix = 0; zix < ids.size( ); zix++ )
         for ( int zp = 0; zp < 2; zp++ )
         {   int z = 2 * ids[zix] + zp;
             double rclock = WallClockTime( );
             basevector r = R[z];

             // Log progress.
             if (LOG_PROGRESS) 
             #pragma omp critical
               dots_pct(reads_processed++, 2 * ids.size());
             
             if (ANNOUNCE) {
               #pragma omp critical
               { cout << "begin read " << z << endl;
               }
             }

             // Start logging.

             ostringstream outx;
             ostream& out = ( DIRECT ? cout : outx );
             if (DIRECT)
               {    out << "\n=========================================================="
                        << "==========================\n\n";    }
             out << "read " << z/2 << ", length = " << R[z].size( ) << endl;

             // Align the read to the genome.

             vec<align_data> adata;
             if ( VALIDATE && READ_LOG != "none" )
               {    double gclock = WallClockTime( );
                 for ( int pass = 1; pass <= 2; pass++ )
                   {    basevector b = r;
                     if ( pass == 2 ) b.ReverseComplement( );
                     vec< triple<int,int,int> > places;
                     for ( int s = 0; s <= b.isize( ) - LG; s++ )
                       {    int n = KmerId( b, LG, s );
                         for ( int z = 0; z < Glocs[n].isize( ); z++ )
                           {    int gid = Glocs[n][z].first;
                             places.push( gid, s, Glocs[n][z].second );    }    }
                     align a;
                     int errors;
                     vec<Bool> checked( places.size( ), False );
                     for ( int q = 0; q < places.isize( ); q++ )
                       {    if ( checked[q] ) continue;
                         int gid = places[q].first;
                         int rpos = places[q].second, gpos = places[q].third;
                         const basevector& g = genome[gid];
     
                         // Test with an initial, cheaper alignment.
     
                         int flank = ( G_initial - LG ) / 2;
                         int b_start = Max( 0, rpos - flank );
                         int b_stop = Min( rpos + LG + flank, b.isize( ) );
                         basevector b0( b, b_start, b_stop - b_start );
                         int bw = G_initial / G_initial_bandwidth_div;
                         int offset = ( rpos - b_start ) - gpos;
                         SmithWatBandedA( b0, g, offset, bw, a, errors, 0, 1, 1 );
                         if ( a.pos1( ) > 0 || a.Pos1( ) < b0.isize( ) ) continue;
                         if ( errors > double( b0.size( ) )
                              * G_initial_error_rate_max )
                           {    continue;    }

                         // Now do the full alignment.
     
                         offset = rpos - gpos;
                         int bandwidth = int(ceil( 
                                                  G_bandwidth_frac * double( r.size( ) ) ));
                         int sub = 1;
                         int ins = 1;
                         int del = 1;
                         int errors2;
                         SmithWatBandedA2<unsigned short>( b, g, offset, 
                                                           bandwidth, a, errors2, (ostream*) 0, sub, ins, del );
                         if ( a.pos1( ) > 0 || a.Pos1( ) < r.isize( ) ) continue;
                         adata.push( a, gid, a.pos2( ), a.Pos2( ), errors2, 
                                     pass == 1 );    

                         // Mark places that are effectively checked by this align.
                         // It might save time to switch to a binary search of 
                         // "places".

                         vec<ho_interval> perfs1, perfs2;
                         a.PerfectIntervals1( b, g, perfs1 );
                         a.PerfectIntervals2( b, g, perfs2 );
                         for ( int w = 0; w < places.isize( ); w++ )
                           {    if ( checked[w] ) continue;
                             if ( places[w].first != gid ) continue;
                             int rpos = places[w].second, gpos = places[w].third;
                             for ( int y = 0; y < perfs1.isize( ); y++ )
                               {    if ( !Subset( ho_interval( rpos, rpos + LG ),
                                                  perfs1[y] ) )
                                   {    continue;    }
                                 if ( !Subset( ho_interval( gpos, gpos + LG ),
                                               perfs2[y] ) )
                                   {    continue;    }
                                 if ( perfs1[y].Start( ) - perfs2[y].Start( )
                                      != rpos - gpos )
                                   {    continue;    }
                                 checked[w] = True;
                                 break;    }    }    }    }

                 UniqueSort(adata);
                 int min_errors = ( adata.empty( ) ? -1 : adata[0].errors );
                 if (PERF_STATS)
                   out << TimeSince(gclock) << " used aligning to genome" << endl;
                 const double max_excess_errors = 0.15;
                 vec<Bool> to_delete_a( adata.size( ), False );
                 for ( int i = 0; i < adata.isize( ); i++ )
                   {    if ( !( adata[i].errors
                                <= double(min_errors) * ( 1.0 + max_excess_errors ) ) )
                       {    to_delete_a[i] = True;    }    }
                 EraseIf( adata, to_delete_a );
                 out << "\n" << adata.size( ) << " genomic placements:\n\n";
                 for ( int i = 0; i < adata.isize( ); i++ )
                   {    int gid = adata[i].gid;
                     out << "[" << i+1 << "] " << gid << "." << adata[i].pos2 << "-" 
                         << adata[i].Pos2 << " " << ( adata[i].fw ? "fw" : "rc" )
                         << ", " << adata[i].errors << " errors\n";
                     basevector b = r;
                     if ( !adata[i].fw ) b.ReverseComplement( );
                     PrintVisualAlignment( False, out, b, genome[gid], adata[i].a );
                   }    }

             // Define alignments of the contigs to the read.  First define local
             // alignments, then associate a global alignment to each local alignment.

             vec< triple<int,int,int> > aligns;
             GetLocalAligns(
                            r, T, Tlocs, L, flank, max_errs, aligns, out, SHOW_ALIGNS, True );
             Sort(aligns);
             vec<align> aligns_a( aligns.size( ) );
             GetGlobalAligns( r, T, aligns, bandwidth_div, aligns_a, 1, 1, 1 );

             // Filter out bad alignments.
                    
             double max_error_rate_remove = 0.37;
             vec<Bool> to_remove( aligns.size( ), False );
             for ( int i = 0; i < aligns_a.isize( ); i++ )
               {    const align& a = aligns_a[i];
                 int u = aligns[i].first;
                 int errs = ActualErrors( T[u], r, a, 1, 1 );
                 double err_rate = double(errs) / double( a.extent1( ) );
                 // PRINT4_TO( out, u, errs, a.extent1( ), err_rate ); // XXXXXXXXXXXX
                 if ( err_rate > max_error_rate_remove )
                   to_remove[i] = True;    }
             EraseIf( aligns, to_remove );
             EraseIf( aligns_a, to_remove );

             // For a given contig, group the alignments.  If we have two local
             // alignments of a read to a given contig, the question we have to
             // answer is whether they extend to the same global alignment, in which
             // case they belong in the same group.

             vec< vec< triple<int,int,int> > > alignsx;
             vec<align> alignsx_a;
             for ( int i = 0; i < aligns.isize( ); i++ )
               {    int j, l, u = aligns[i].first;
                 vec< triple<int,int,int> > xxx;
                 for ( j = i; j < aligns.isize( ); j++ )
                   {    if ( aligns[j].first != aligns[i].first ) break;
                     xxx.push_back( aligns[j] );    }
                 equiv_rel e( xxx.size( ) );
                 for ( int k = 0; k < xxx.isize( ); k++ )
                   {    int rpos = xxx[k].second;
                     int offset = xxx[k].third;
                     int predicted_overlap = IntervalOverlap( 
                                                             0, T[u].isize( ), offset, offset + r.isize( ) );
                     int bandwidth = predicted_overlap / bandwidth_div;
                     align a;
                     int errors;
                     SmithWatBandedA(T[u], r, offset, bandwidth, a, errors, 0, 1, 1);
                     vec<ho_interval> perfs1, perfs2;
                     a.PerfectIntervals1( T[u], r, perfs1 );
                     a.PerfectIntervals2( T[u], r, perfs2 );
                     vec<Bool> to_remove( perfs1.size( ), False );
                     for ( int l = 0; l < perfs1.isize( ); l++ )
                       if ( perfs1[l].Length( ) < L ) to_remove[l] = True;
                     EraseIf( perfs1, to_remove );
                     EraseIf( perfs2, to_remove );
                     for ( int l = 0; l < xxx.isize( ); l++ )
                       {    if ( l == k ) continue;
                         int rpos = xxx[l].second;
                         int upos = xxx[l].third + rpos;
                         for ( int r = 0; r < perfs1.isize( ); r++ )
                           {    if ( !Member( perfs1[r], upos ) ) continue;
                             if ( !Member( perfs2[r], rpos ) ) continue;
                             if ( upos - perfs1[r].Start( ) 
                                  != rpos - perfs2[r].Start( ) )
                               {    continue;    }
                             e.Join( k, l );    }    }    }
                 vec<int> reps;
                 e.OrbitRepsAlt(reps);
                 for ( int l = 0; l < reps.isize( ); l++ )
                   {    vec<int> o;
                     e.Orbit( reps[l], o );
                     vec< triple<int,int,int> > x;
                     for ( int r = 0; r < o.isize( ); r++ )
                       x.push_back( xxx[ o[r] ] );
                     sort( x.begin( ), x.end( ), cmp_third );
                     alignsx.push_back(x);
                     alignsx_a.push_back( aligns_a[ i + o[0] ] );    }
                 i = j - 1;    }

             // Create indices.

             int NN = alignsx.size( );
             vec<int> tot(NN);
             for ( int i = 0; i < NN; i++ )
               tot[i] = alignsx[i][0].first;

             // Print alignments.

             for ( int i = 0; i < NN; i++ )
             {    const align& a = alignsx_a[i];
                  int t = tot[i];
                  out << "contig " << tot[i] << ", pos1 = " << a.pos1( ) << "-" 
                      << a.Pos1( ) << " of " << T[t].size( ) << ", pos2 = "
                      << a.pos2( ) << "-" << a.Pos2( ) << " of "
                      << R[z].size( ) << "\n";    }

             // Look for potential gap closers.

             int npatches = 0;
             for ( int i1 = 0; i1 < NN; i1++ )
               {    for ( int i2 = 0; i2 < NN; i2++ )
                   {    int t1 = tot[i1], t2 = tot[i2];
                     int s = to_super[t1];
                     if ( to_super[t2] != s ) continue;
                     if ( to_super_pos[t2] != to_super_pos[t1] + 1 ) continue;
                     const align &a1 = alignsx_a[i1], &a2 = alignsx_a[i2];
                     if ( a1.Pos1( ) != T[t1].isize( ) ) continue;
                     if ( a2.pos1( ) != 0 || a1.pos2( ) != 0 ) continue;
                     if ( a2.Pos2( ) != R[z].isize( ) ) continue;

                     /*
                     if ( a1.extent1( ) < 200 ) continue;
                     if ( a2.extent1( ) < 200 ) continue;
                     */

                     // Build a gap patcher.  Note a few issues here:
                     // (1) The farther we have to reach back on a read to find
                     // an L-mer match, the worse the read.  We might want to
                     // discriminate against the patches that reach far back.
                     // (2) If there are multiple alignments we take the first one.
                     // This can't be optimal.

                     out << "building a gap patcher" << endl;
                     vec< triple<int,int,int> > av1 = alignsx[i1];
                     vec< triple<int,int,int> > av2 = alignsx[i2];
                     sort( av1.begin( ), av1.end( ), cmp_second );
                     sort( av2.begin( ), av2.end( ), cmp_second );

                     // Compute the 'effective' length of the alignments.

                     const double max_err = 0.05;
                     basevector rrc(r);
                     rrc.ReverseComplement( );
                     align a1rc(a1);
                     a1rc.ReverseThis( T[t1].size( ), r.size( ) );
                     int best_len1 = BestGoodFirst( a1rc, Trc[t1], rrc, max_err );
                     int best_len2 = BestGoodFirst( a2, T[t2], r, max_err );
                     int best_len = Min( best_len1, best_len2 );

                     // Choose anchors that are a bit back from the ends
                     // of the unipaths.

                     const int leeway = 20;
                     vec<ho_interval> perf_tig1, perf_read1;
                     a1.PerfectIntervals1( T[t1], r, perf_tig1 );
                     a1.PerfectIntervals2( T[t1], r, perf_read1 );
                     vec<ho_interval> perf_tig2, perf_read2;
                     a2.PerfectIntervals1( T[t2], r, perf_tig2 );
                     a2.PerfectIntervals2( T[t2], r, perf_read2 );
                     int n1 = T[t1].size( ), j1, j2;
                     for ( j1 = perf_tig1.isize( ) - 1; j1 >= 0; j1-- )
                     {    if ( perf_tig1[j1].Length( ) < L ) continue;
                          if ( n1 - perf_tig1[j1].Start( ) >= leeway + L )
                              break;    }
                     if ( j1 < 0 )
                     {    out << "Couldn't find left anchor\n";
                          continue;    }
                     int tpos1 = perf_tig1[j1].Start( );
                     int rpos1 = perf_read1[j1].Start( );
                     for ( int x = 0; x < perf_tig1[j1].Length( ) - L; x++ )
                     {    if ( n1 - tpos1 == leeway + L ) break;
                          rpos1++, tpos1++;    }
                     for ( j2 = 0; j2 < perf_tig2.isize( ); j2++ )
                     {    if ( perf_tig2[j2].Length( ) < L ) continue;
                          if ( perf_tig2[j2].Stop( ) >= leeway ) break;    }
                     if ( j2 == perf_tig2.isize( ) )
                     {    out << "Couldn't find right anchor\n";
                          continue;    }
                     int tpos2 = perf_tig2[j2].Stop( );
                     int rpos2 = perf_read2[j2].Stop( );
                     for ( int x = 0; x < perf_tig2[j2].Length( ); x++ )
                     {    if ( tpos2 == leeway ) break;
                          rpos2--, tpos2--;    }

                     // Build the patch.

                     if ( rpos1 < rpos2 + L )
                     {   basevector rtrim( r, rpos1, rpos2 + L - rpos1 );
                         int expected_gap = a2.pos2( ) - a1.Pos2( );
                         int read_id = z/2;
                         #pragma omp critical
                         {    if (PATCH_STATS)
                              {    out << "\n";
                                   PRINT6_TO( out, t1, t2, read_id, rtrim, 
                                        tpos1, tpos2 );
                                   int lerrs = ActualErrors( T[t1], r, a1, 1, 1 );
                                   int rerrs = ActualErrors( T[t2], r, a2, 1, 1 );
                                   int lmerrs = ActualErrors( T[t1], r, a1, 1, 0 );
                                   int rmerrs = ActualErrors( T[t2], r, a2, 1, 0 );
                                   double err_rate = double( lerrs + rerrs )
                                        / double( a1.extent2( ) + a2.extent2( ) );
                                   double merr_rate = double( lmerrs + rmerrs )
                                        / double( a1.extent2( ) + a2.extent2( ) );
                                   int min_over = Min( a1.extent2(), a2.extent2() );
                                   PRINT6_TO( out, min_over, expected_gap, 
                                        a1.pos1( ), a2.Pos1( ), 
                                        err_rate, merr_rate );
                                   PRINT2_TO( out, rpos1, rpos2 );    
                                   PRINT2_TO( out, best_len1, best_len2 );
                                   out << "\nleft alignment:\n";
                                   PrintVisualAlignment( False, out, T[t1], r, a1 );
                                   out << "right alignment:\n";
                                   PrintVisualAlignment( 
                                        False, out, T[t2], r, a2 );    }
                              patchers.push( t1, t2, rtrim, z, tpos1, 
                                   tpos2 + L, best_len, expected_gap );    }
                         npatches++;    }    }    }

             // Add to log.

             if ( ( npatches > 0 && READ_LOG != "none" ) || READ_LOG == "all" )
               {    if (PERF_STATS)
                   {    out << "\ntime used on this read = " << TimeSince(rclock)
                            << "\n";    }
                 report[z] = outx.str( );    }    }

       // Sort and filter patchers.  For the filtering, we first computed "best_len"
       // above.  It is defined as follows.  For each alignment (of the read to
       // the left unipath and the right unipath), we compute the maximum number
       // of bases on the contig end that are in alignment with substitution rate 
       // <= 5%.  Then we take the minimum (over the two alignments).  This is
       // best_len.  Then for a given gap, we reverse order the patches by best_len,
       // and take the top 10 (or top half, if less).  For these, we compute the 
       // median expected_gap.  Then we go through all the patches for the gap, and 
       // discard those that are not in the top group and are not within 10% of the 
       // median.

       UniqueSort(patchers);
       vec<Bool> patchers_to_delete( patchers.size( ), False );
       for ( int i = 0; i < patchers.isize( ); i++ )
       {    int j;
            for ( j = i + 1; j < patchers.isize( ); j++ )
            {    if ( patchers[j].t1 != patchers[i].t1 ) break;
                 if ( patchers[j].t2 != patchers[i].t2 ) break;    }
            if ( j - i < 2 )
            {    i = j - 1;
                 continue;    }
            vec< pair<int,int> > bleg;
            for ( int k = i; k < j; k++ )
                 bleg.push( patchers[k].best_len, patchers[k].expected_gap );
            ReverseSort(bleg);
            int keep = Min( 10, bleg.isize( )/2 );
            int worst = bleg[ keep - 1 ].first;
            vec<int> ex;
            for ( int k = 0; k < keep; k++ )
                 ex.push_back( bleg[k].second );
            Sort(ex);
            int mid = ex[ ex.isize( ) / 2 ];
            // PRINT3( patchers[i].t1, patchers[i].t2, mid ); // XXXXXXXXXXXXXXXXXXX
            double delta = 0.1;
            int low = int( floor( (1.0 - delta) * double(mid) ) );
            int high = int( ceil( (1.0 + delta) * double(mid) ) );
            vec<int> used;
            for ( int k = i; k < j; k++ )
            {    if ( Member( used, patchers[k].rid ) ) 
                      patchers_to_delete[k] = True;
                 else used.push_back( patchers[k].rid );
                 if ( patchers[k].best_len < worst
                      && ( patchers[k].expected_gap < low 
                        || patchers[k].expected_gap > high ) )
                 {    patchers_to_delete[k] = True;    }    }
            i = j - 1;    }
       EraseIf( patchers, patchers_to_delete );

       // Print reports.

       int rcount = 0;
       for ( size_t i = 0; i < R.size( ); i++ )
         {    if ( report[i].size( ) == 0 ) continue;
           if ( !DIRECT )
             {    cout << "\n=========================================================="
                       << "==========================\n\n";    }
           cout << "<" << ++rcount << "> " << report[i];    }
       flush(cout);
     if ( !PATCH ) return 0;
     
       // ---- Writing patches to disk
       cout << Tag() << "Writing patchers to '" << fn_patchers << "'." << endl;
       BinaryWriter::writeFile(fn_patchers.c_str(), patchers);
       cout << Tag() << "Wrote " << patchers.size() << " patches." << endl;
     }



     
     // ---- Reading patches from disk

     cout << Tag() << "Reading patchers from '" << fn_patchers << "'." << endl;
     vec<GapPatcher> patchers;
     BinaryReader::readFile(fn_patchers.c_str(), & patchers);
     const size_t n_patchers = patchers.size();
     cout << Tag() << "Read " << n_patchers << " patches." << endl;
     UniqueSort(patchers);


     // ---- Subdivide patchers by gaps

     const int nb_rad_SW = 10;  // the radius of the local Smith-Waterman
     const unsigned n_patchers_min = 5;
     const unsigned sz_padding_min = 2 * nb_rad_SW + L;
     vec< vec<GapPatcher> > patchers_gap;
     patchers_gap_collect(T, patchers, n_patchers_min, sz_padding_min,
                          &patchers_gap);
     const size_t n_gaps = patchers_gap.size();


     // ---- Choose, for each gap, the first-guess patcher 
     
     const unsigned type = 0; // 0:shortest size  1:median size  2:median gap
     const unsigned sz_patcher_min = 2 * L;
     vec<unsigned> ips_best;
     patchers_first_guesses(T, patchers_gap, type, sz_patcher_min, 
                            & ips_best);


     // ---- Sort gaps according to their predicted performance
     //      with omp use schedule(dynamic, 1) to run largest jobs first 

     vec<unsigned> is_gaps;
     indexes_gaps_for_parallel(patchers_gap, ips_best,
                               & is_gaps);


     // ---- Declare some stuff and optimize the guess patchers


     cout << "\n" << Tag() << ": patching (100 dots to follow)" << endl;

     vec<String> yreports(n_gaps);
     vec<bool> patched(ntigs, false);
     vec<GapPatcher0> all_paths(ntigs);

     uint n_gaps_done = 0;
     vec<double> timers_sw(3, 0);

     #pragma omp parallel for schedule(dynamic, 1)
     for (unsigned j_gap = 0; j_gap < n_gaps; j_gap++) {

       const unsigned i_gap = is_gaps[j_gap];
       const vec<GapPatcher> & patchers = patchers_gap[i_gap];
       const unsigned n_patchers = patchers.size();

       const unsigned ip_best = ips_best[i_gap];
       if (ip_best >= n_patchers) {  // invalid best
         if (PATCH_VERBOSITY >= 1) 
           #pragma omp critical
           cout << Date() << ": i_gap= " << i_gap << ": invalid best patch." << endl; 
       }
       else {
         const size_t    ibv1 = patchers[0].t1;
         const size_t    ibv2 = patchers[0].t2;
         
         GapPatcher0 & p_opt = all_paths[ibv1];
         vec<double> timers(3, 0);
         patcher_optimal(T, patchers, ip_best, L, nb_rad_SW, sz_padding_min, i_gap, 
                         & p_opt, PATCH_VERBOSITY, & timers);
         patched[ibv1] = true;

         #pragma omp critical
         { 
           timers_sw[0] += timers[0];
           timers_sw[1] += timers[1];
           timers_sw[2] += timers[2];
         }


         ostringstream out_tmp;

         if (VALIDATE_PATCHES) {
           const BaseVec & bv1  = T[ibv1];
           const BaseVec & bv2  = T[ibv2];
           BaseVec c = Cat(basevector(bv1, 0, p_opt.upos1), 
                           p_opt.r, 
                           basevector(bv2, p_opt.upos2, bv2.isize() - p_opt.upos2));
           vec<placementx> totals = FindGenomicPlacements(c, LG, genome, Glocs);
           #pragma omp critical
           {
             for (int l = 0; l < totals.isize(); l++) {
               out_tmp << "perfectly matches " << totals[l].g << "."
                       << totals[l].pos << "-"
                       << totals[l].pos + c.isize() << "\n";
             }
             if ( totals.empty( ) ) {
               Ofstream( out, "slobber.fasta" );                         
               c.Print( out, "patch_for_gap_" + ToString(j_gap) );
               out_tmp << AllOfOutput1( "QueryLookupTable K=12 MM=12 MC=0.15 "
                                        "SEQS=slobber.fasta L=" + data_dir + "/genome.lookup "
                                        "VISUAL=True NH=True QUIET=True" ); 
             }    

           }
         }
         
         yreports[i_gap] = out_tmp.str();

         // -- Log progress.

         #pragma omp critical
         { if (LOG_PROGRESS) dots_pct(n_gaps_done++, n_gaps); }
       }
     }



     cout << endl;
     cout << Tag() << ": times for smith-waterman:" << endl 
          << setw(10) << timers_sw[0] << " secs total" << endl
          << setw(10) << timers_sw[1] << " secs short" << endl
          << setw(10) << timers_sw[2] << " secs long" << endl;

     int yc = 1;
     for (int i = 0; i < yreports.isize(); i++) {
       if (yreports[i].size() > 0) {
         cout << "\n[" << yc++ << "] " << yreports[i];
       }
     }

     // Save patches, after minimizing them so that the patches themselves are as 
     // small as possible.

     if (WRITE)
     {    cout << "\n" << Tag( ) << ": saving patches" << endl;
          vec<assembly_edit> edits;
          for ( int s = 0; s < scaffolds.isize( ); s++ ) 
          {    const superb& S = scaffolds[s];
               for ( int i = 0; i < S.Ngaps( ); i++ )
               {    int m1 = S.Tig(i), m2 = S.Tig(i+1);
                    if ( patched[m1] )
                    {    int start1 = all_paths[m1].upos1;
                         int stop2 = all_paths[m1].upos2;
                         basevector r = all_paths[m1].r;
                         while( start1 < (int) tigsa[m1].size( )
                              && r.size( ) > 0
                              && tigsa[m1][start1] == as_base( r[0] ) )
                         {    start1++;
                              r.SetToSubOf( r, 1, r.isize( ) - 1 );    }
                         while( stop2 > 0 && r.size( ) > 0
                              && tigsa[m2][stop2-1] == as_base( r[ r.isize()-1 ] ) )
                         {    stop2--;
                              r.resize( r.isize( ) - 1 );    }
                         vec<basevector> reps;
                         reps.push_back(r);
                         edits.push( assembly_edit::GAP_CLOSER, m1, 
                              start1, m2, stop2, reps );    }    }    }
          BinaryWriter::writeFile( outputFile, edits );    }
     cout << Sum(patched) << " gaps patched" << endl;
     cout << endl << Tag() << "Done, time used = " << TimeSince(clock) << endl;    }
