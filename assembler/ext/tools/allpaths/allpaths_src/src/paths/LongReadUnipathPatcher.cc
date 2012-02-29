//////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// LongReadUnipathPatcher.  Align long reads to the unipath graph.
// Work in progress.

// To do:
// 1. Don't write reads that are subsumed by unipaths.
// 3. Speed up.
// 4. Get rid of MAX_EDGES and MAX_PATHS arguments.
// 5. Why do we report multiple genomic placements for reads, at the same position?
// 6. In some cases we report success, but in fact we've found an unnecessary gap.
// 7. Don't use the trimmed and aligned reads.  Trimming must be de novo.  Also
//    aligning puts all the reads in the 'forward' orientation.
// 8. Handle the case where the graph has cycles.
// 9. Use read groups (from splitting at smart bells) better and don't
//    discard the tiny reads associated with them.
// 10. Improve handling of gap consensus
//     (a) definition of individual gap sequences
//     (b) method for forming consensus
//     (c) may need to allow for alternatives.
// 11. Do better job testing correctness of gap sequences.
// 12. Are we losing half the reads?  Consensus piles are not fw-rc symmetric.
// 13. Should test to see what happens when some reads are randomly reverse
//     complemented.  It is likely to expose bugs.
// 14. Method of joining can't be right because it will always accept a join of
//     min_overlap size even if the true sequence is a gap.
// 15. Association of a full alignment to each local alignment by GetGlobalAligns: 
//     very inefficient and also redundant with next step.
// 16. More accurately estimate error profile.
// 17. Don't find all paths.  Instead find an acyclic covering graph.
// 18. Now finding cycles, didn't use to.
// 19. DistX uses SmithWatFree, which is inefficient: should use a bounded 
//     Smith-Waterman with both ends anchored.
// 20. Switched SmithWatBandedA2 parameters to sub = ins = del = 1 because I'm not
//     sure it's working.
// 21. Note that templatization of SmithWatBandedA2 is completely stupid.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "Equiv.h"
#include "FastIfstream.h"
#include "FetchReads.h"
#include "MainTools.h"
#include "ParseSet.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/GetNexts.h"
#include "paths/LongReadTools.h"
#include "paths/LongReadUnipathPatcherUtils.h"
#include "paths/UnibaseUtils.h"
#include "random/Random.h"

double gclock = WallClockTime( );

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

// An edgepath is a path through a directed graph.  It specifies vertex ids and
// edge ids in the graph.

class edgepath {

     public:

     edgepath( ) { }
     edgepath( const vec<int>& verts, const vec<int>& edges )
          : verts(verts), edges(edges)
     {    ForceAssertEq( verts.size( ), edges.size( ) + 1 );    }

     vec<int> verts, edges;

};

class uniseq {

     public:

     static vecbasevector* unibasesp;

     uniseq( ) { }
     uniseq( const vec<int>& unis, const vec<int>& overs ) : unis(unis), overs(overs)
     {    ForceAssertEq( unis.size( ), overs.size( ) + 1 );    }

     String ToStr( ) const
     {    String s = ToString( unis[0] );
          for ( int i = 0; i < overs.isize( ); i++ )
          {    s += " --[" + ToString(overs[i]) + "]--> " + ToString(unis[i+1]);    }
          return s;    }

     basevector Bases( ) const
     {    basevector b = (*unibasesp)[ unis[0] ];
          for ( int j = 0; j < overs.isize( ); j++ )
          {    b.resize( b.isize( ) - overs[j] );
               b = Cat( b, (*unibasesp)[ unis[j+1] ] );    }
          return b;    }

     friend Bool operator==( const uniseq& u1, const uniseq& u2 )
     {    return u1.unis == u2.unis && u1.overs == u2.overs;    }
     friend Bool operator!=( const uniseq& u1, const uniseq& u2 )
     {    return !( u1 == u2 );    }
     friend Bool operator<( const uniseq& u1, const uniseq& u2 )
     {    if ( u1.unis < u2.unis ) return True;
          if ( u1.unis > u2.unis ) return False;
          return u1.overs < u2.overs;    }
     friend Bool operator>( const uniseq& u1, const uniseq& u2 )
     {    if ( u1.unis > u2.unis ) return True;
          if ( u1.unis < u2.unis ) return False;
          return u1.overs > u2.overs;    }

     friend uniseq Join( const uniseq& u1, const uniseq& u2, const int overlap )
     {    vec<int> unis = u1.unis;
          unis.append( u2.unis );
          vec<int> overs = u1.overs;
          overs.push_back(overlap);
          overs.append( u2.overs );
          return uniseq( unis, overs );    }

     vec<int> unis;
     vec<int> overs;

     void Reverse( const vec<int>& to_rc )
     {    unis.ReverseMe( );
          overs.ReverseMe( );
          for ( int i = 0; i < unis.isize( ); i++ )
               unis[i] = to_rc[ unis[i] ];    }

};

vecbasevector* uniseq::unibasesp;

void SmithWatFreeSym( const basevector& b1, const basevector& b2, align& a )
{    alignment al;
     int best_loc;
     if ( b1.size( ) <= b2.size( ) )
     {    SmithWatFree( b1, b2, best_loc, al, false, false, 1, 1 );
          a.UnpackFrom(al);    }
     else
     {    SmithWatFree( b2, b1, best_loc, al, false, false, 1, 1 );
          a.UnpackFrom(al);
          a.Flip( );    }    }



// A gap patcher is defined by left and right uniseqs u1 and u2, and a basevector r,
// whose left end aligns to u1 starting at upos1 and whose right end aligns to u2
// ending at upos2, where the positions are in terms of the basevectors associated
// to u1 and u2.  The idea is that the ends of r would have been trimmed back to
// places where it matches u1 and u2 along an L-mer.
//
//                                 upos1             upos2
//                                 *                 *
//    -----------------u1---------------    ---------------u2----------------
//                                 -------r-----------

class gap_patcher {

     public:

     gap_patcher( ) { }
     gap_patcher( const uniseq& u1, const uniseq& u2, const basevector& r,
          const int rid, const int upos1, const int upos2 )
          : u1(u1), u2(u2), r(r), rid(rid), upos1(upos1), upos2(upos2) { }

     friend Bool operator==( const gap_patcher& p1, const gap_patcher& p2 )
     {    return p1.u1 == p2.u1 && p1.u2 == p2.u2 && p1.r == p2.r
               && p1.rid == p2.rid && p1.upos1 == p2.upos1;    }

     friend Bool operator<( const gap_patcher& p1, const gap_patcher& p2 )
     {    if ( p1.u1 < p2.u1 ) return True;
          if ( p1.u1 > p2.u1 ) return False;
          if ( p1.u2 < p2.u2 ) return True;
          if ( p1.u2 > p2.u2 ) return False;
          if ( p1.r < p2.r ) return True;
          if ( p1.r > p2.r ) return False;
          if ( p1.rid < p2.rid ) return True;
          if ( p1.rid > p2.rid ) return False;
          if ( p1.upos1 < p2.upos1 ) return True;
          if ( p1.upos1 > p2.upos1 ) return False;
          return p1.upos2 < p2.upos2;    }

     uniseq u1, u2;
     basevector r;
     int rid;
     int upos1, upos2;

     void Reverse( const vec<int>& to_rc )
     {    basevector b1 = u1.Bases( ), b2 = u2.Bases( );
          u1.Reverse(to_rc), u2.Reverse(to_rc);
          swap( u1, u2 );
          r.ReverseComplement( );
          int upos1_new = b2.isize( ) - upos2, upos2_new = b1.isize( ) - upos1;
          upos1 = upos1_new, upos2 = upos2_new;    }

};




int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     // Data set.
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);

     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
          "Number of threads to use (use all available processors if set to 0)");

     CommandArgument_String_OrDefault(IN_HEAD, "all_reads");
     CommandArgument_String_OrDefault(OUT_HEAD, "all_reads.patched");
     CommandArgument_String_OrDefault(READS, "pacbio_reads_orig");

     CommandArgument_String_OrDefault(PACBIO_RUNS, "");
     CommandArgument_String_OrDefault_Doc(IDS, "all", 
          "list of reads in ParseIntSet format or \"all\"");
     CommandArgument_Int_OrDefault_Doc(MIN_READ_SIZE, 150,
          "ignore reads shorter than this");
     // Algorithmic heuristics.
     CommandArgument_Bool_OrDefault(NEUTER, False);
     CommandArgument_Int_OrDefault_Doc(MAX_EDGES, 1000,
          "give up if more than this many edges in graph");
     CommandArgument_Int_OrDefault_Doc(MAX_PATHS, 10000,
          "give up if more than this many paths in graph");
     CommandArgument_Bool_OrDefault_Doc(PATCH_ONLY, False, 
          "'True' skips reads that can not be used for patching");
     CommandArgument_Bool_OrDefault_Doc(TERMINAL_ONLY, False, 
          "only align to unipaths that have a dead end");
     // Evaluation.
     CommandArgument_Bool_OrDefault_Doc(VALIDATE, False,
	  "assess results versus reference genome");
     CommandArgument_String_OrDefault(GENOME_FASTB_PATH, "");
     CommandArgument_Bool_OrDefault(COMPUTE_TRUE_READ_LENGTH, False);
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
     // Saving.
     CommandArgument_Bool_OrDefault_Doc(WRITE, True, 
          "Save new unibases");
     EndCommandArguments;

     // Thread control
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define directories.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String unibasefile = run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K);
     String unipathfile = run_dir + "/" + IN_HEAD + ".unibases.k" + ToString(K);
     cout << Date( ) << ": " << run_dir << endl;

     // Check for null operation.

     String KS = ToString(K);
     if (NEUTER)
     {    vecbasevector U(unibasefile);
          U.WriteAll( run_dir + "/" + OUT_HEAD + ".unibases.k"+ KS );   
          return 0;    }

     // Load data.

     vecbasevector R0, genome;
     if (PACBIO_RUNS != "") {
       vec<int> runs;
       ParseIntSet( PACBIO_RUNS, runs );
       for ( int i = 0; i < runs.isize( ); i++ ) {
	 vecbasevector x;
	 FetchReads( x, 0, "/seq/pacbio_results/"
		     "smrtanalysis-1.0.1/common/jobs/017/0" + ToString( runs[i] )
		     + "/data/filtered_subreads.fa" );
	   R0.Append(x);
       }
     } else {
       R0.ReadAll( run_dir + "/" + READS + ".fastb");
     }

     if (VALIDATE) 
     {    String g = ( GENOME_FASTB_PATH == "" ? data_dir + "/genome.fastb"
               : GENOME_FASTB_PATH );
          genome.ReadAll( g );    }

     // Define heuristic constants for validation.

     const int LG = 12;
     const int validation_twiddle = 10;

     // Define heuristic constants for de novo algorithm.

     const int L = 12;
     const int flank = 19;
     const int max_errs = 10;
     double max_error_rate = 0.26;
     const double bandwidth_div = 5;

     // Set up.

     vecbasevector U( unibasefile );
     uniseq::unibasesp = &U;
     vec<int> to_rc;
     UnibaseInvolution( U, to_rc, K );
     vec< vec<int> > nexts;
     GetNexts( K, U, nexts );
     if (TERMINAL_ONLY)
     {    for ( size_t i = 0; i < U.size( ); i++ )
          {    if ( nexts[i].nonempty( ) && nexts[ to_rc[i] ].nonempty( ) )
                    U[i].resize(0);    }    }
     ForceAssert( READ_LOG == "all" || READ_LOG == "hard" || READ_LOG == "none" );

     // Load PathsToLocs output file.

     vec<String> ptl;
     vec< triple<int,int,int> > ptl_loc;

     // Validation needs PathToLocs (not available here)
     // if (VALIDATE)
     // {    fast_ifstream in( run_dir + "/all_reads.unipaths.k" + ToString(K)
     //           + ".PathsToLocs.out" );
     //      String line;
     //      while(1)
     //      {    getline( in, line );
     //           if ( in.fail( ) ) break;
     //           if ( line.Contains( "path ", 0 ) )
     //           {    ptl.push_back(line);    
     //                int tig = line.Between( "> ", "." ).Int( );
     //                int start = line.Between( ".", "-" ).Int( );
     //                int stop = line.After( ">" ).Between( "-", " " ).Int( );
     //                ptl_loc.push( tig, start, stop );    }    }    }

     // Report reads stats and remove very short ones.

     vecbasevector R;
     longlong total_length = 0, sub_total = 0;
     for ( size_t i = 0; i < R0.size( ); i++ )
     {    total_length += R0[i].size( );
          if ( R0[i].isize( ) >= MIN_READ_SIZE ) 
          {    R.push_back( R0[i] );
               sub_total += R0[i].size( );    }    }
     cout << "have " << ToStringAddCommas( R0.size( ) ) 
          << " reads of mean length " << total_length / R0.size( ) << "\n";
     cout << "using " << ToStringAddCommas( R.size( ) ) 
          << " reads of length >= " << MIN_READ_SIZE;
     if (VALIDATE)
     {    longlong gbases = 0;
          for ( size_t i = 0; i < genome.size( ); i++ )
               gbases += genome[i].size( );
          cout << ", and providing " << setprecision(3) 
               << double(sub_total)/double(gbases) << "x coverage";    }
     cout << endl;
     vec<int> ids;
     if ( IDS == "all" )
     {    for ( size_t j = 0; j < R.size( ); j++ )
               ids.push_back(j);    }
     else
     {    ParseIntSet( IDS, ids, True, 0, R.size( ) );
          for ( int j = 0; j < ids.isize( ); j++ )
               ForceAssertLe( ids[j], (int) R.size( ) );    }

     // Hash the genome.

     vec< vec< pair<int,int> > > Glocs;
     if (VALIDATE) BuildGLocs( LG, genome, Glocs );

     // Hash the unipaths.

     vec< vec< pair<int,int> > > Ulocs( IPow( 4, L ) );
     for ( size_t i = 0; i < U.size( ); i++ )
     {    for ( int j = 0; j <= U[i].isize( ) - L; j++ )
          {    int n = KmerId( U[i], L, j );
               Ulocs[n].push( i, j );    }    }

     // Set up global stats.

     int simple_good = 0, found_alls = 0;

     // Describe error profile.

     /*
     vec<int> MGG( 3, 0 );
     int nbases = 0;
     const int read_sample = 100;
     #pragma omp parallel for
     for ( int i = 0; i < read_sample; i++ )
     {    int z = randomx( ) % R.size( );
          const basevector& r = R[z];
          vec< triple<int,int,int> > aligns;
          GetLocalAligns( r, U, Ulocs, L, flank, max_errs, aligns, cout );
          Sort(aligns);
          vec<align> aligns_a( aligns.size( ) );
          GetGlobalAligns( r, U, aligns, bandwidth_div, aligns_a, 1, 1, 1 );
          for ( int j = 0; j < aligns.isize( ); j++ )
          {    const align& a = aligns_a[j];
               int u = aligns[j].first;
               if ( a.pos2( ) > 0 || a.Pos2( ) < r.isize( ) ) continue;
               vec<int> mgg = a.MutationsGap1Gap2( U[u], r );
               #pragma omp critical
               {    for ( int l = 0; l < 3; l++ )
                         MGG[l] += mgg[l];
                    nbases += r.size( );    }
               break;    }    }
     double sub_frac = double(MGG[0]) / double(nbases);
     double ins_frac = double(MGG[1]) / double(nbases);
     double del_frac = double(MGG[2]) / double(nbases);
     cout << "\nsubstitutions = " << setiosflags(ios::fixed) << setprecision(1) 
          << 100.0 * sub_frac << "%, insertions = " << 100.0 * ins_frac << "%, "
          << "deletions = " << 100.0 * del_frac << "%" 
          << resetiosflags(ios::fixed) << endl;
     PRINT(nbases);
     */
     double sub_frac = 0.05, ins_frac = 0.05, del_frac = 0.05; // not really used

     // Make an initial pass through the reads to find local alignments to the
     // unibases.

     cout << Date( ) <<  ": finding local alignments (100 dots)" << endl;
     int reads_processed = 0;
     vec< vec< triple<int,int,int> > > ALIGNS( ids.size( ) );
     int dots_printed = 0;
     #pragma omp parallel for
     for ( size_t zi = 0; zi < ids.size( ); zi++ )
     {    int z = ids[zi];
          double rclock = WallClockTime( );
          basevector r = R[z];
          #pragma omp critical
          {    reads_processed++;    
               if (LOG_PROGRESS)
               {    int newdots 
                         = ( 100 * reads_processed ) / ids.size( )
                         - ( 100 * (reads_processed-1) ) / ids.size( );
                    if ( newdots > 0 )
                    {    for ( int j = 0; j < newdots; j++ )
                         {    cout << ".";
                              dots_printed++;
                              if ( dots_printed % 50 == 0 ) cout << "\n";
                              else if ( dots_printed % 10 == 0 ) cout << " ";    }
                         flush(cout);    }    }    }
          ostringstream out;
          GetLocalAligns( r, U, Ulocs, L, flank, max_errs, ALIGNS[zi], out, 
                          SHOW_ALIGNS, PATCH_ONLY );    }

     // Look for reads that need to be split at inversions.  These cases correspond
     // to undetected or degenerate smartbells.  We specifically exclude the case 
     // where a palindromic inversion is present in the genome.

     cout << Date( ) << ": splitting reads at undetected smartbells" << endl;
     vecbasevector R2;
     vec< vec< triple<int,int,int> > > ALIGNS2;
     vec<String> ORIGIN;
     for ( size_t zi = 0; zi < ids.size( ); zi++ )
     {    int z = ids[zi];
          vec< pair<int,int> > breakpoints;
          for ( int m = 0; m < ALIGNS[zi].isize( ) - 1; m++ )
          {    int u = ALIGNS[zi][m].first;
               if ( to_rc[u] == ALIGNS[zi][m+1].first )
               {    if ( U[u].Overlap( U[ to_rc[u] ], K-1 ) ) continue;
                    int p1 = ALIGNS[zi][m].second;
                    int p2 = ALIGNS[zi][m+1].second;
                    breakpoints.push( p1+L, p2 );

                    /*
                    cout << "\nSEE EVIDENCE OF INVERSION\n";
                    for ( int r = m; r <= m+1; r++ )
                    {    int i = ALIGNS[zi][r].first, p = ALIGNS[zi][r].second;
                         int j = ALIGNS[zi][r].second + ALIGNS[zi][r].third;
                         cout << "i = " << i << ", j = " << j << ", p = " << p
                              << ", j-p = " << j-p << ", U[i].size = " 
                              << U[i].size( ) << "\n";    }    
                    */

                    }    }

          for ( int j = 0; j < breakpoints.isize( ); j++ )
          {    int start = ( j == 0 ? 0 : breakpoints[j-1].second );
               int stop = breakpoints[j].first;
               basevector r( R[z], start, stop - start );
               R2.push_back_reserve(r);
               ORIGIN.push_back( ToString(z) + ":" + ToString(j+1) );
               vec< triple<int,int,int> > v;
               for ( int l = 0; l < ALIGNS[zi].isize( ); l++ )
               {    if ( ALIGNS[zi][l].second < start ) continue;
                    if ( ALIGNS[zi][l].second + L > stop ) continue;
                    v.push( ALIGNS[zi][l].first, ALIGNS[zi][l].second - start,
                         ALIGNS[zi][l].third + start );    }
               ALIGNS2.push_back(v);    }
          int start = ( breakpoints.empty( ) ? 0 : breakpoints.back( ).second );
          int stop = R[z].isize( );
          basevector r( R[z], start, stop - start );
          R2.push_back_reserve(r);
          if ( breakpoints.empty( ) ) ORIGIN.push_back( ToString(z) );
          else 
          {    ORIGIN.push_back( 
                    ToString(z) + ":" + ToString( breakpoints.size( ) + 1 ) );    }
          vec< triple<int,int,int> > v;
          for ( int l = 0; l < ALIGNS[zi].isize( ); l++ )
          {    if ( ALIGNS[zi][l].second < start ) continue;
               if ( ALIGNS[zi][l].second + L > stop ) continue;
               v.push( ALIGNS[zi][l].first, ALIGNS[zi][l].second - start,
                    ALIGNS[zi][l].third + start );    }
          ALIGNS2.push_back(v);    }


     // Sort indices based on ALIGNS2[i].size() * R2[i].size()
     vec<size_t> iR_sorted(R2.size());
     {
       const size_t nr = R2.size();
       vec< pair<size_t, size_t> >  pairs(nr);
       for (unsigned ir = 0; ir != nr; ir++) {
         pairs[ir].first = ALIGNS2[ir].size() * R2[ir].size();
         pairs[ir].second = ir;
       }
       sort(pairs.rbegin(), pairs.rend());
       for (size_t ir = 0; ir != nr; ir++) {
         iR_sorted[ir] = pairs[ir].second;
       }
     } 
     



     // Go through the reads.
     vec<gap_patcher> patchers;
     vec< triple<int,int,int> > newtigs;
     vec<String> report( R2.size( ) );
     vecbasevector all_paths;

     reads_processed = 0;
     dots_printed = 0;
     cout << Date( ) << ": main pass through reads (100 dots)" << endl;
     #pragma omp parallel for schedule(dynamic,1)
     for ( size_t jR2 = 0; jR2 < R2.size( ); jR2++ )
     {    double rclock = WallClockTime( );
          const size_t iR2 = iR_sorted[jR2];
          basevector r = R2[iR2];

          // Log progress.

          #pragma omp critical
          {    reads_processed++;    
               if (LOG_PROGRESS)
               {    int newdots 
                         = ( 100 * reads_processed ) / R2.size( )
                         - ( 100 * (reads_processed-1) ) / R2.size( );
                    if ( newdots > 0 )
                    {    for ( int j = 0; j < newdots; j++ )
                         {    cout << ".";
                              dots_printed++;
                              if ( dots_printed % 50 == 0 ) cout << "\n";
                              else if ( dots_printed % 10 == 0 ) cout << " ";    }
                         flush(cout);    }    }    }
          if (ANNOUNCE)
          {
               #pragma omp critical
               {    cout << "begin read " << ORIGIN[iR2] << endl;    }    }

          // Start logging.

          ostringstream outx;
          ostream& out = ( DIRECT ? cout : outx );
          if (DIRECT)
          {    out << "\n=========================================================="
                    << "==========================\n\n";    }
          out << "read " << ORIGIN[iR2] << ", length = " << R2[iR2].size( ) << endl;

          // Align the read to the genome.

          vec<align_data> adata;
          if ( VALIDATE && READ_LOG != "none" ) {
	    AlignToGenome( LG, COMPUTE_TRUE_READ_LENGTH,
			   r, genome, Glocs, adata, out );
	    if (PERF_STATS)
	      out << TimeSince(gclock) << " used aligning to genome" << endl;
	  }

          // Define alignments of the unipaths to the read.  First define local
          // alignments (done above), then associate a global alignment to each 
          // local alignment.

          Sort(ALIGNS2[iR2]);
          vec<align> aligns_a( ALIGNS2[iR2].size( ) );
          GetGlobalAligns( r, U, ALIGNS2[iR2], bandwidth_div, aligns_a, sub_frac,
               ins_frac, del_frac );

          // Filter out bad alignments.
                    
          double max_error_rate_remove = 0.37;
          vec<Bool> to_remove( ALIGNS2[iR2].size( ), False );
          for ( int i = 0; i < aligns_a.isize( ); i++ )
          {    const align& a = aligns_a[i];
               int u = ALIGNS2[iR2][i].first;
               int errs = ActualErrors( U[u], r, a, 1, 1 );
               double err_rate = double(errs) / double( a.extent1( ) );
               // PRINT4_TO( out, u, errs, a.extent1( ), err_rate ); // XXXXXXXXXXXX
               if ( err_rate > max_error_rate_remove )
                    to_remove[i] = True;    }
          EraseIf( ALIGNS2[iR2], to_remove );
          EraseIf( aligns_a, to_remove );

          /*
          for ( int i = 0; i < aligns_a.isize( ); i++ )
          {    const align& a = aligns_a[i];
               int u = ALIGNS2[iR2][i].first;
               int errs = ActualErrors( U[u], r, a, 1, 1 );
               double err_rate = double(errs) / double( a.extent1( ) );
               PRINT2_TO( out, u, err_rate );    }
          */

          // Find the unipaths that are touched by the read.

          vec<int> us;
          for ( int i = 0; i < ALIGNS2[iR2].isize( ); i++ )
               us.push_back( ALIGNS2[iR2][i].first );
          UniqueSort(us);
          int nu = us.size( );

          // For a given unipath, group the alignments.  If we have two local
          // alignments of a read to a given unipath, the question we have to
          // answer is whether they extend to the same global alignment, in which
          // case they belong in the same group.

          vec< vec< triple<int,int,int> > > alignsx;
          vec<align> alignsx_a;
          for ( int i = 0; i < ALIGNS2[iR2].isize( ); i++ )
          {    int j, l, u = ALIGNS2[iR2][i].first;
               vec< triple<int,int,int> > xxx;
               for ( j = i; j < ALIGNS2[iR2].isize( ); j++ )
               {    if ( ALIGNS2[iR2][j].first != ALIGNS2[iR2][i].first ) break;
                    xxx.push_back( ALIGNS2[iR2][j] );    }
               // sort( xxx.begin( ), xxx.end( ), cmp_second );
               equiv_rel e( xxx.size( ) );
               for ( int k = 0; k < xxx.isize( ); k++ )
               {    int rpos = xxx[k].second;
                    int offset = xxx[k].third;
                    int predicted_overlap = IntervalOverlap( 
                         0, U[u].isize( ), offset, offset + r.isize( ) );
                    int bandwidth = predicted_overlap / bandwidth_div;
                    align a;
                    int errors;
                    SmithWatBandedA(U[u], r, offset, bandwidth, a, errors, 0, 1, 1);
                    vec<ho_interval> perfs1, perfs2;
                    a.PerfectIntervals1( U[u], r, perfs1 );
                    a.PerfectIntervals2( U[u], r, perfs2 );
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
          vec<int> tou(NN), touid(NN);
          for ( int i = 0; i < NN; i++ )
          {    tou[i] = alignsx[i][0].first;
               touid[i] = BinPosition( us, tou[i] );    }

          // Find the implied overlaps.

          vec< triple<int,int,int> > u1u2overlap;
          for ( int i1 = 0; i1 < ALIGNS2[iR2].isize( ); i1++ )
          {    for ( int i2 = 0; i2 < ALIGNS2[iR2].isize( ); i2++ )
               {    if ( i1 == i2 ) continue;
                    if ( ALIGNS2[iR2][i2].second != ALIGNS2[iR2][i1].second ) 
                         continue;
                    int u1 = ALIGNS2[iR2][i1].first, u2 = ALIGNS2[iR2][i2].first;

                    // Compute the start position of u2 relative to u1.
          
                    int offset = ALIGNS2[iR2][i1].third - ALIGNS2[iR2][i2].third;
                    if ( offset < 0 ) continue;

                    // Compute the implied overlap and test for illegal overlap.

                    int overlap = IntervalOverlap(
                         -offset, -offset + U[u1].isize( ), 0, U[u2].isize( ) );
                    if ( !PerfectMatch(U[u1], U[u2], offset, 0, overlap) ) continue;

                    // Save overlap.
     
                    u1u2overlap.push( u1, u2, overlap );    }    }
          UniqueSort(u1u2overlap);

          // Identify "left" and "right" alignment groups, according to whether
          // the unipath appears to extend the read to the left (or the right).
          // These are approximate statements that could probably be made more
          // accurate by aligning rather than inferring start and stop points, but
          // they would still be approximate.

          vec<Bool> left( NN, False ), right( NN, False );
          for ( int i = 0; i < NN; i++ )
          {    for ( int j = 0; j < alignsx[i].isize( ); j++ )
               {    int u = alignsx[i][j].first;
                    int offset = alignsx[i][j].third;
                    if ( j == 0 && offset >= 0 ) left[i] = True;
                    if ( j == alignsx[i].isize( ) - 1 
                         && offset + r.isize( ) <= U[u].isize( ) ) 
                    {    right[i] = True;    }    }    }

          // Build a graph G whose vertices are alignment groups, and whose edges
          // are implied overlaps, as defined by u1u2overlap.

          vec< vec<int> > fromx(NN), tox(NN), from_edge_objx(NN), to_edge_objx(NN);
          vec<int> edgesx;
          for ( int i1 = 0; i1 < NN; i1++ )
          {    int u1 = alignsx[i1][0].first;
               for ( int i2 = 0; i2 < NN; i2++ )
               {    if ( i1 == i2 ) continue;
                    int u2 = alignsx[i2][0].first;
                    vec<int> overs;
                    for ( int r1 = 0; r1 < alignsx[i1].isize( ); r1++ )
                    {    for ( int r2 = 0; r2 < alignsx[i2].isize( ); r2++ )
                         {    if ( alignsx[i1][r1].second != alignsx[i2][r2].second )
                                   continue;
                              int offset 
                                   = alignsx[i1][r1].third - alignsx[i2][r2].third;
                              if ( offset < 0 ) continue;
                              int overlap = IntervalOverlap( -offset, 
                                   -offset + U[u1].isize( ), 0, U[u2].isize( ) );
                              if ( !Member( u1u2overlap,
                                   make_triple( u1, u2, overlap ) ) )
                              {    continue;    }
                              overs.push_back(overlap);    }    }
                    UniqueSort(overs);
                    for ( int r = 0; r < overs.isize( ); r++ )
                    {    fromx[i1].push_back(i2);
                         tox[i2].push_back(i1);
                         from_edge_objx[i1].push_back( edgesx.size( ) );
                         to_edge_objx[i2].push_back( edgesx.size( ) );
                         edgesx.push_back( overs[r] );    }    }    }
          for ( int i = 0; i < NN; i++ )
          {    SortSync( fromx[i], from_edge_objx[i] );
               SortSync( tox[i], to_edge_objx[i] );    }
          digraphE<int> G( fromx, tox, edgesx, to_edge_objx, from_edge_objx );

          // If a sink and a source overlap by K-1, connect them.

          vec<int> sources, sinks;
          G.Sources(sources), G.Sinks(sinks);
          for ( int i1 = 0; i1 < sinks.isize( ); i1++ )
          {    for ( int i2 = 0; i2 < sources.isize( ); i2++ )
               {    int u1 = tou[ sinks[i1]] , u2 = tou[ sources[i2] ];
                    if ( U[u1].Overlap( U[u2], K-1 ) )
                         G.AddEdge( sinks[i1], sources[i2], K-1 );    }    }

          // There are two passes.  In the first pass, we eliminate all 'incomplete'
          // edges.  If the first pass fails, we do a second pass in which all edges
          // are used.

          digraphE<int> G0(G);
          vec<int> to_delete0;
          for ( int i = 0; i < G0.EdgeObjectCount( ); i++ )
               if ( G0.EdgeObject(i) < K - 1 ) to_delete0.push_back(i);
          G0.DeleteEdges(to_delete0);
          G0.RemoveDeadEdgeObjects( );
          for ( int pass = 1; pass <= 2; pass++ )
          {    out << "STARTING PASS " << pass << endl;
               digraphE<int>& GG = ( pass == 1 ? G0 : G );

               // Eliminate edges that are not needed.

               vec<Bool> to_delete( GG.EdgeObjectCount( ), False );
               for ( int v = 0; v < GG.N( ); v++ )
               {    for ( int j = 0; j < GG.From(v).isize( ); j++ )
                    {    int w = GG.From(v)[j];
                         int e = GG.EdgeObjectIndexByIndexFrom( v, j );
                         if ( to_delete[e] ) continue;
                         if ( GG.Source(v) && GG.Sink(w) ) continue;
     
                         // We will try to eliminate e.  First check to see if all
                         // compositions of e with edges to the right of it are 
                         // present.
     
                         Bool OK = True;
                         for ( int k = 0; k < GG.From(w).isize( ); k++ )
                         {    int x = GG.From(w)[k];
                              int f = GG.EdgeObjectIndexByIndexFrom( w, k );
                              if ( to_delete[f] ) continue;

                              // Compute the composition of e with f and see if it's
                              // present.
     
                              int o = GG.EdgeObject(e) + GG.EdgeObject(f) 
                                   - U[ tou[w] ].isize( );
                              Bool found = False;
                              for ( int l = 0; l < GG.From(v).isize( ); l++ )
                              {    if ( GG.From(v)[l] != x ) continue;
                                   int h = GG.EdgeObjectIndexByIndexFrom( v, l );
                                   if ( to_delete[h] ) continue;
                                   if ( GG.EdgeObjectByIndexFrom( v, l ) == o )
                                   {    found = True;
                                        break;    }    }
                              if ( !found )
                              {    OK = False;
                                   break;    }    }
                         if ( !OK ) continue;

                         // Now check to see if compositions of e with edges to the
                         // left of it are present.
     
                         for ( int k = 0; k < GG.To(v).isize( ); k++ )
                         {    int x = GG.To(v)[k];
                              int f = GG.EdgeObjectIndexByIndexTo( v, k );
                              if ( to_delete[f] ) continue;
     
                              // Compute the composition of f with e and see if it's
                              // present.

                              int o = GG.EdgeObject(e) + GG.EdgeObject(f) 
                                   - U[ tou[v] ].isize( );
                              Bool found = False;
                              for ( int l = 0; l < GG.From(x).isize( ); l++ )
                              {    if ( GG.From(x)[l] != w ) continue;
                                   int h = GG.EdgeObjectIndexByIndexFrom( x, l );
                                   if ( to_delete[h] ) continue;
                                   if ( GG.EdgeObjectByIndexFrom( x, l ) == o )
                                   {    found = True;
                                        break;    }    }
                              if ( !found )
                              {    OK = False;
                                   break;    }    }
                         if (OK) to_delete[e] = True;    }    }
               vec<int> to_deletex;
               for ( int e = 0; e < to_delete.isize( ); e++ )
                    if ( to_delete[e] ) to_deletex.push_back(e);
               GG.DeleteEdges(to_deletex);
     
               // Print the graph.

               digraphE<int> GGred(GG);
               GGred.RemoveDeadEdgeObjects( );
               out << "\nGRAPH HAS " << GGred.EdgeObjectCount( ) << " EDGES\n";
               if ( GGred.EdgeObjectCount( ) > MAX_EDGES )
               {    out << "MAX_EDGES exceeded" << "\n";
                    out << "Giving up." << "\n";
                    GG.Clear( );    }
               int last_u1 = -1;
               vec<int> whichu( GG.N( ) );
               {    vec<int> counter( us.size( ), 0 );
                    for ( int x = 0; x < GG.N( ); x++ )
                         whichu[x] = ++counter[ touid[x] ];    }
               for ( int x1 = 0; x1 < GG.N( ); x1++ )
               {    int u1 = tou[x1];
                    if ( u1 != last_u1 ) 
                    {    if ( GG.From(x1).nonempty( ) ) out << "\n";
                         last_u1 = u1;    }
                    for ( int j = 0; j < GG.From(x1).isize( ); j++ )
                    {    int x2 = GG.From(x1)[j];
                         int u2 = tou[x2];
                         int overlap = GG.EdgeObjectByIndexFrom( x1, j );
                         if ( left[x1] ) out << "L-";
                         out << "u1 = " << u1 << "[" << whichu[x1] << "], u2 = "
                              << u2 << "[" << whichu[x2] << "]";
                         if ( right[x2] ) out << "-R";
                         out << ", overlap = " << overlap << "\n";    }    }
               for ( int x = 0; x < GG.N( ); x++ )
               {    if ( GG.From(x).empty( ) && GG.To(x).empty( ) )
                    {    out << "\n";
                         if ( left[x] ) out << "L-";
                         out << tou[x] << "[" << whichu[x] << "]";
                         if ( right[x] ) out << "-R";
                         out << "\n";    }    }

               // Find all paths through GG.  Note that there are computational 
               // performance issues with this if reads are sufficiently long
               // and polymorphism is sufficiently great.
               //
               // Warning: this preliminary version does not work if the graph has
               // cycles.
     
               vec<basevector> bpaths;
               vec<uniseq> uniseqs;
               vec< vec<int> > edge_paths;
               vec<int> to_left, to_right;
               Bool fail = False;
               if ( !GG.Acyclic( ) )
                    out << "Graph has cycles, won't find all paths." << endl;
               else
               {    vec<int> sources, sinks;
                    GG.Sources(sources), GG.Sinks(sinks);
                    for ( int i1 = 0; i1 < sources.isize( ); i1++ )
                    {    for ( int i2 = 0; i2 < sinks.isize( ); i2++ )
                         {    if ( !left[ sources[i1] ] && !right[ sinks[i2] ] )
                                   continue;
                              vec< vec<int> > edge_paths_this;
                              GG.EdgePaths( sources[i1], sinks[i2], 
                                   edge_paths_this, -1, MAX_PATHS );
                              edge_paths.append(edge_paths_this);    
                              if ( edge_paths.isize( ) > MAX_PATHS )
                              {    out << "\nMAX_PATHS exceeded" << "\n";
                                   out << "Giving up." << "\n";
                                   GG.Clear( );    
                                   edge_paths.clear( );
                                   fail = True;
                                   break;    }    }
                         if (fail) break;    }
                    if (fail)
                    {    if ( READ_LOG != "none" ) report[iR2] = outx.str( );
                         break;    }
                    out << "\n";
                    if (PERF_STATS) PRINT_TO( out, edge_paths.size( ) );
                    if (PERF_STATS) out << Date( ) << ": building bpaths" << endl;
                    out << "PATHS\n";
                    GG.ToLeft(to_left), GG.ToRight(to_right);
                    for ( int i = 0; i < edge_paths.isize( ); i++ )
                    {    const vec<int>& p = edge_paths[i];
                         basevector b = U[ tou[ to_left[ p[0] ] ] ];
                         for ( int j = 0; j < p.isize( ); j++ )
                         {    basevector c = U[ tou[ to_right[ p[j] ] ] ];
                              c.SetToSubOf( c, GG.EdgeObject(p[j]), 
                                   c.isize( ) - GG.EdgeObject(p[j]) );
                              b = Cat( b, c );    }
                         bpaths.push_back(b);    }
                    if (PERF_STATS) out << Date( ) << ": sorting bpaths" << endl;
                    UniqueSortSync( bpaths, edge_paths );
                    if (PERF_STATS) out << Date( ) << ": sort complete" << endl;
                    if (PERF_STATS) PRINT_TO( out, edge_paths.size( ) );
                    for ( int i = 0; i < edge_paths.isize( ); i++ )
                    {    out << "\npath " << i+1 << ":";
                         const vec<int>& p = edge_paths[i];
                         out << " " << tou[ to_left[ p[0] ] ];
                         for ( int j = 0; j < p.isize( ); j++ )
                         {    out << " --[" << GG.EdgeObject(p[j]) << "]--> "
                                   << tou[ to_right[ p[j] ] ];    }
                         out << "\n";    
                         vec<int> unis, overs;
                         unis.push_back( tou[ to_left[ p[0] ] ] );
                         for ( int j = 0; j < p.isize( ); j++ )
                         {    overs.push_back( GG.EdgeObject(p[j]) );
                              unis.push_back( tou[ to_right[ p[j] ] ] );     }
                         uniseqs.push( unis, overs );

                         // Find perfect matches to the genome.

                         if ( VALIDATE && READ_LOG != "none" )
                         {    PrintGenomicPlacements( 
                                   out, bpaths[i], L, genome, Glocs );    }    }    }

               // Add missing unipaths.  After this bpaths is in bijective 
               // correspondence with edge_paths + vert_paths.

               vec<int> vert_paths;
               {    vec<Bool> used( us.size( ), False );
                    for ( int i = 0; i < edge_paths.isize( ); i++ )
                    {    const vec<int>& p = edge_paths[i];
                         for ( int j = 0; j < p.isize( ); j++ )
                         {    used[ touid[ to_left[ p[j] ] ] ] = True;
                              used[ touid[ to_right[ p[j] ] ] ] = True;    }    }
                    int path_count = edge_paths.size( );
                    for ( int i = 0; i < GG.N( ); i++ )
                    {    int uid = touid[i];
                         if ( used[uid] ) continue;
                         vert_paths.push_back(i);
                         used[uid] = True;
                         bpaths.push_back( U[ tou[i] ] );    
                         out << "\npath " << ++path_count << ": " << tou[i] << "\n";
                         if ( VALIDATE && READ_LOG != "none" )
                         {    PrintGenomicPlacements( out, bpaths.back( ), L, 
                                   genome, Glocs );    }    
                         vec<int> unis, overs;
                         unis.push_back( tou[i] );
                         uniseqs.push( unis, overs );    }    }

               // Define LL and RR.

               vec<Bool> LL( bpaths.size( ) ), RR( bpaths.size( ) );
               for ( int i = 0; i < bpaths.isize( ); i++ )
               {    if ( i < edge_paths.isize( ) )
                    {    LL[i] = left[ to_left[ edge_paths[i].front( ) ] ];
                         RR[i] = right[ to_right[ edge_paths[i].back( ) ] ];    }
                    else
                    {    LL[i] = left[ vert_paths[ i - edge_paths.isize( ) ] ];
                         RR[i] = right[ vert_paths[ i - edge_paths.isize( ) ] ];
                             }    }

               // Make edgepaths.  This combines edge_paths + vert_paths.

               vec<edgepath> edgepaths;
               for ( int i = 0; i < edge_paths.isize( ); i++ )
               {    vec<int> verts, edges = edge_paths[i];
                    verts.push_back( to_left[ edges.front( ) ] );
                    for ( int j = 0; j < edges.isize( ); j++ )
                         verts.push_back( to_right[ edges[j] ] );
                    edgepaths.push( verts, edges );    }
               for ( int i = 0; i < vert_paths.isize( ); i++ )
               {    vec<int> verts, edges;
                    verts.push_back( vert_paths[i] );
                    edgepaths.push( verts, edges );    }

               // At this point the data structures in use are edgepaths, bpaths, 
               // uniseqs, L, and R, which are all in parallel.  We try to make 
               // some joins.

               const int min_overlap = 10;
               Bool have_LR = False;
               int nb = bpaths.size( );
               for ( int i = 0; i < nb; i++ )
                    if ( LL[i] && RR[i] ) have_LR = True;
               if ( !have_LR )
               {    for ( int i1 = 0; i1 < nb; i1++ )
                    for ( int i2 = 0; i2 < nb; i2++ )
                    {    if ( !LL[i1] || !RR[i2] ) continue;
                         int v1 = edgepaths[i1].verts.back( );
                         int v2 = edgepaths[i2].verts.front( );
                         out << "\nconsidering overlaps between " << tou[v1] << "[" 
                              << whichu[v1] << "] and "
                              << tou[v2] << "[" << whichu[v2] << "]\n";
                         int u1 = tou[v1], u2 = tou[v2];
                         const align &a1 = alignsx_a[v1], &a2 = alignsx_a[v2];

                         /*
                         out << "first alignment:\n";
                         PrintVisualAlignment( False, out, U[u1], r, a1 );
                         out << "second alignment:\n";
                         PrintVisualAlignment( False, out, U[u2], r, a2 );
                         */

                         int gap = a2.pos2( ) - a1.Pos2( );
                         if ( gap >= 0 )
                              out << "predict gap of " << gap << " bases\n";
                         else out << "predict overlap of " << -gap << " bases\n";
                         for ( int o = min_overlap; o <= K-2; o++ )
                         {    if ( Overlap( bpaths[i1], bpaths[i2], o ) )
                              {    out << "see overlap of " << o << " bases\n";
                                   basevector x = bpaths[i1];
                                   x.resize( x.isize( ) - o );
                                   x = Cat( x, bpaths[i2] );
                                   bpaths.push_back(x);
                                   uniseqs.push_back( 
                                        Join( uniseqs[i1], uniseqs[i2], o ) );
                                   LL.push_back(True); 
                                   RR.push_back(True);    }    }    }    }

               // Create L-R paths.  These are either one L-R bpath (in which case
               // second is -1), or an L-/-R pair, if there are no L-R bpaths.
               // Note that in the L-/-R case we concatenate the basevectors, which
               // doesn't make a whole lot of sense.

               vec< triple<int,int,basevector> > lrpaths;
               for ( int pass = 1; pass <= 2; pass++ ) 
               {    for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
                    {    if ( pass == 2 && lrpaths.nonempty( ) ) break;
                         if ( pass == 1 )
                         {    if ( LL[i1] && RR[i1] )
                              {    lrpaths.push( i1, -1, bpaths[i1] );
                                   if ( lrpaths.isize( ) > MAX_PATHS ) 
                                        break;    }
                              continue;    }
                         if ( !LL[i1] ) continue;
                         for ( int i2 = 0; i2 < bpaths.isize( ); i2++ )
                         {    if ( !LL[i2] && RR[i2] )
                              {    lrpaths.push( i1, i2, Cat( bpaths[i1], 
                                        bpaths[i2] ) );    
                                   if ( lrpaths.isize( ) > MAX_PATHS ) 
                                        break;    }
                              if ( lrpaths.isize( ) > MAX_PATHS ) 
                                   break;    }    }    }
               if ( lrpaths.isize( ) > MAX_PATHS )
               {    out << "\nMAX_PATHS exceeded" << "\n";
                    out << "Giving up." << "\n";
                    fail = True;
                    if ( READ_LOG != "none" ) report[iR2] = outx.str( );
                    break;    }

               // Smith-Waterman the read against each L-R path, report results.
               // Select those within 15% of best.

               out << "\nFINAL ALIGNS\n" << endl;
               vec< triple<int,basevector,int> > lrpath_aligns;
               vec<int> left_trim, right_trim;
               for ( int i = 0; i < lrpaths.isize( ); i++ )
               {    if ( r.size( ) > lrpaths[i].third.size( ) ) continue;
                    int best_loc;
                    alignment a;
                    // run time could be reduced by speeding this up:
                    int errs = SmithWatFree( r, lrpaths[i].third, best_loc, a,
                         false, false, 1, 1 );
                    basevector b( lrpaths[i].third, a.pos2( ), a.extent2( ) );
                    lrpath_aligns.push( errs, b, i );
                    left_trim.push_back( a.pos2( ) );
                    right_trim.push_back( 
                         lrpaths[i].third.isize( ) - a.Pos2( ) );    }
               SortSync( lrpath_aligns, left_trim, right_trim );
               int min_errors = 1000000000;
               for ( int i = 0; i < lrpath_aligns.isize( ); i++ )
                    min_errors = Min( min_errors, lrpath_aligns[i].first );
               double error_rate = double(min_errors) / double( r.isize( ) );
               const double max_excess_errors = 0.15;
               int count = 0;
               Bool sg = False;
               vec<Bool> found( adata.size( ), False );
               for ( int i = 0; i < lrpath_aligns.isize( ); i++ )
               {    int j;
                    for ( j = i + 1; j < lrpath_aligns.isize( ); j++ )
                    {    if ( lrpath_aligns[j].second != lrpath_aligns[i].second ) 
                              break;    }
                    if ( lrpath_aligns[i].first 
                         > double(min_errors) * ( 1.0 + max_excess_errors ) )
                    {    i = j - 1;
                         continue;    }
                    int lrid = lrpath_aligns[i].third;
                    int l1 = lrpaths[lrid].first, l2 = lrpaths[lrid].second;
                    String bid = ( l2 < 0 ? ToString(l1+1)
                         : ToString(l1+1) + "+" + ToString(l2+1) );
                    lrpath_aligns[i].second.Print( out, "[" + ToString(++count) 
                         + "] path " + bid + ", "
                         + ToString( lrpath_aligns[i].first ) + " errors" );
                    if ( l2 >= 0 )
                    {    basevector b1 = bpaths[l1], b2 = bpaths[l2];
                         Bool bad = False;
                         if ( left_trim[i] >= b1.isize( ) ) bad = True;
                         if ( right_trim[i] >= b2.isize( ) ) bad = True;
                         int v1 = edgepaths[l1].verts.back( );
                         int v2 = edgepaths[l2].verts.front( );

                         // Require that there is a dead end on at least one side.

                         if ( nexts[ tou[v1] ].nonempty( ) 
                              || nexts[ to_rc[ tou[v2] ] ].nonempty( ) )
                         {    bad = True;    }

                         if ( !bad )
                         {    
                              // Build a gap patcher.  

                              out << "building a gap patcher" << endl;
                              vec< triple<int,int,int> > av1 = alignsx[v1];
                              vec< triple<int,int,int> > av2 = alignsx[v2];
                              sort( av1.begin( ), av1.end( ), cmp_second );
                              sort( av2.begin( ), av2.end( ), cmp_second );

                              // Choose anchors that are a bit back from the ends
                              // of the unipaths.

                              int m1 = av1.isize( ) - 1, m2 = 0;
                              const int leeway = 20;
                              while ( m1 > 0 && U[ tou[v1] ].isize( ) 
                                   - av1[m1].second - av1[m1].third - L < leeway ) 
                              {    m1--;    }
                              while( m2 < av2.isize( ) - 1 
                                   && av2[m2].second + av2[m2].third < leeway )
                              {    m2++;    }
                              PRINT4_TO(out, m1, m2, av1.size(), av2.size()); // XXX

                              // Build the patch.

                              int rpos1 = av1[m1].second, rpos2 = av2[m2].second;
                              if ( rpos1 < rpos2 + L )
                              {    basevector rtrim( r, rpos1, rpos2 + L - rpos1 );
                                   int upos1 = av1[m1].third + rpos1
                                        + bpaths[l1].isize( ) 
                                        - U[ tou[v1] ].isize( );
                                   int upos2 = av2[m2].third + rpos2 + L;
                                   #pragma omp critical
                                   {    patchers.push( uniseqs[l1], uniseqs[l2],
                                             rtrim, iR2, 
                                             upos1, upos2 );    }    }    }    }

                    if (VALIDATE)
                    {    int lid = lrpath_aligns[i].third;
                         int bid1 = lrpaths[lid].first, bid2 = lrpaths[lid].second;
                         if ( bid2 < 0 )
                         {    basevector b = lrpath_aligns[i].second;
                              vec<placementx> places = 
                                   FindGenomicPlacements( b, L, genome, Glocs );
                              for ( int j = 0; j < places.isize( ); j++ )
                              {    const placementx& p = places[j];
                                   out << ( p.fw ? "fw" : "rc" ) << " matches "
                                        << "genome " << p.g << "." << p.pos << "-" 
                                        << p.pos + b.isize( ) << " perfectly\n";    
                                   #pragma omp critical
                                   {    newtigs.push( 
                                             p.g, p.pos, p.pos + b.isize( ) );    }
                                   for ( int q = 0; q < adata.isize( ); q++ )
                                   {    const align_data& aq = adata[q];
                                        if ( aq.gid != p.g ) continue;
                                        if ( aq.fw != p.fw ) continue;
                                        if ( Abs( aq.pos2 - p.pos ) 
                                             > validation_twiddle ) 
                                        {    continue;    }
                                        if ( Abs( aq.Pos2 - p.pos - b.isize( ) )
                                             > validation_twiddle )
                                        {    continue;    }
                                        out << "matches genomic placement "
                                             << q+1 << "\n";    
                                        found[q] = True;    }    }    }
                         else
                         {    basevector b1 = bpaths[bid1], b2 = bpaths[bid2];
                              vec<placementx> places1, places2;
                              if ( left_trim[i] < b1.isize( ) 
                                   && right_trim[i] < b2.isize( ) )
                              {    b1.SetToSubOf( b1, left_trim[i],
                                        b1.isize( ) - left_trim[i] );
                                   b2.resize( b2.isize( ) - right_trim[i] );
                                   places1 = FindGenomicPlacements( 
                                        b1, L, genome, Glocs );
                                   places2 = FindGenomicPlacements( 
                                        b2, L, genome, Glocs );    }
                              for ( int j1 = 0; j1 < places1.isize( ); j1++ )
                              for ( int j2 = 0; j2 < places2.isize( ); j2++ )
                              {    const placementx& p1 = places1[j1];
                                   const placementx& p2 = places2[j2];
                                   if ( p1.g != p2.g || p1.fw != p2.fw ) continue;
                                   int gap, start, stop;
                                   if (p1.fw) 
                                   {    gap = p2.pos - p1.pos - b1.isize( );
                                        start = p1.pos;    }
                                   else 
                                   {    gap = p1.pos - p2.pos - b2.isize( );
                                        start = p2.pos;    }
                                   if ( gap < 0 ) continue;
                                   stop = start + b1.isize( ) + b2.isize( ) + gap;
                                   out << ( p1.fw ? "fw" : "rc" ) << " matches "
                                        << "genome " << p1.g << "." << start << "-" 
                                        << stop
                                        << " perfectly with gap " << gap << "\n";
                                   for ( int q = 0; q < adata.isize( ); q++ )
                                   {    const align_data& aq = adata[q];
                                        if ( aq.gid != p1.g ) continue;
                                        if ( aq.fw != p1.fw ) continue;
                                        if ( Abs( aq.pos2 - start ) 
                                             > validation_twiddle ) 
                                        {    continue;    }
                                        if ( Abs( aq.Pos2 - stop )
                                             > validation_twiddle )
                                        {    continue;    }
                                        out << "matches genomic placement "
                                             << q+1 << "\n";    
                                        found[q] = True;    }    }    }    }
                    i = j - 1;    }
               Bool found_all = False;
               if ( VALIDATE && READ_LOG != "none" )
               {    for ( int w = 0; w <= adata.isize( ); w++ )
                    {    if ( w == adata.isize( ) 
                              || adata[w].errors > adata[0].errors )
                         {    found_all = True;
                              break;    }
                         if ( !found[w] ) break;    }
                    sg = ( count == 1 && found_all );    }

               // Dump report.

               if ( READ_LOG != "none" ) report[iR2] = outx.str( );

               // Are we done?

               if ( error_rate <= max_error_rate ) 
               {    if (VALIDATE)
                    {    if ( found_all && !sg )
                         {    
                              #pragma omp critical
                              {    found_alls++;    }
                              out << "\nfound all\n";
                              if ( READ_LOG != "none" ) report[iR2] = outx.str( );    }
                         if (sg) 
                         {    
                              #pragma omp critical
                              {    simple_good++;    }
                              if ( READ_LOG == "all" ) out << "\nsimple good\n";
                              else report[iR2].clear( );    }    }
                    break;    }    }

          // Dump report.

          if ( READ_LOG != "none" && ( !simple_good || READ_LOG == "all" ) )
          {    if (PERF_STATS)
               {    out << "\ntime used on this read = " << TimeSince(rclock) 
                         << "\n";    }
               report[iR2] = outx.str( );    }
          if (ANNOUNCE)
          {
               #pragma omp critical
               {    cout << "end read " << iR2 << endl;    }    }    }

     // Print reports.

     cout << Date( ) << ": main pass complete" << endl;
     int rcount = 0;
     for ( size_t i = 0; i < R2.size( ); i++ )
     {    if ( report[i].size( ) == 0 ) continue;
          if ( !DIRECT )
          {    cout << "\n=========================================================="
                    << "==========================\n\n";    }
          cout << "<" << ++rcount << "> " << report[i];    }

     // Canonicalize orientation of patchers.

     for ( int i = 0; i < patchers.isize( ); i++ )
     {    gap_patcher x = patchers[i];
          x.Reverse(to_rc);
          if ( x < patchers[i] ) patchers[i] = x;    }

     // Derive consensus from patchers.

     UniqueSort(patchers);
     vec<int> pstarts, pstops;
     for ( int i = 0; i < patchers.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < patchers.isize( ); j++ )
          {    if ( patchers[j].u1 != patchers[i].u1 ) break;
               if ( patchers[j].u2 != patchers[i].u2 ) break;    }
          pstarts.push_back(i), pstops.push_back(j);
          i = j - 1;    }
     vec<String> yreports( pstarts.size( ) );
     #pragma omp parallel for
     for ( int q = 0; q < pstarts.isize( ); q++ )
     {    int start = pstarts[q], stop = pstops[q];
          const int min_patches = 2;
          if ( stop - start < min_patches ) continue;
          ostringstream xout;
          PRINT2_TO( xout, q, pstarts.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          xout << "finding consensus:\n";
          xout << "left = " << patchers[start].u1.ToStr( ) << "\n";
          xout << "right = " << patchers[start].u2.ToStr( ) << "\n";
          basevector b1 = patchers[start].u1.Bases( );
          basevector b2 = patchers[start].u2.Bases( );
          PRINT2_TO( xout, b1.size( ), b2.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXX

          // Define the pack of patches.

          vec<GapPatcher0> p;
          vec<int> prid;
          for ( int j = start; j < stop; j++ )
          {    p.push( patchers[j].r, patchers[j].upos1, patchers[j].upos2 );
               prid.push_back( patchers[j].rid );    }
          SortSync( p, prid );
          int np = p.size( );
          PRINT_TO( xout, np ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          xout << "\npatches:\n";
          for ( int j = 0; j < p.isize( ); j++ )
          {    const int pflank = 30;
               int lstart = Max( 0, p[j].upos1 - pflank ), lstop = p[j].upos1;
               int rstart = p[j].upos2; 
               int rstop = Min( b2.isize( ), p[j].upos2 + pflank );
               basevector left( b1, lstart, lstop - lstart );
               basevector right( b2, rstart, rstop - rstart );
               basevector rp = Cat( left, p[j].r, right );
               xout << ">[" << j + 1 << "] " << ORIGIN[ prid[j] ] << " "
                    << p[j].upos1 << ", " << p[j].upos2 << "\n"
                    << rp.ToString( ) << endl;    }

          // Find the shortest patch p1, however for now don't use patches that
          // are shorter than 2*L, which would break the code.

          int min_patch = 1000000000, best_id = -1;
          for ( int i = 0; i < np; i++ )
          {    if ( p[i].r.isize( ) < min_patch && p[i].r.isize( ) >= 2*L )
               {    min_patch = p[i].r.size( );
                    best_id = i;    }    }
          if ( best_id < 0 ) continue;

          xout << "shortest patch is " << best_id+1 << "\n"; // XXXXXXXXXXXXXXXXXXXX
          PRINT2_TO( xout, p[best_id].upos1, p[best_id].upos2 ); // XXXXXXXXXXXXXXXX

          GapPatcher0 X = patcher_optimal_original(p, best_id, L, b1, b2);


          basevector c = Cat( basevector( b1, 0, X.upos1 ), X.r, 
                              basevector( b2, X.upos2, b2.isize( ) - X.upos2 ) );
          #pragma omp critical
          {    all_paths.push_back_reserve( c );    }
          c.Print( xout, "consensus" );
          if (VALIDATE)
          {    vec<placementx> totals = FindGenomicPlacements( c, LG, genome, Glocs );
               #pragma omp critical
               {    for ( int l = 0; l < totals.isize( ); l++ )
                    {    xout << "perfectly matches " << totals[l].g << "."
                              << totals[l].pos << "-"
                              << totals[l].pos + c.isize( ) << "\n";
                         newtigs.push( totals[l].g, totals[l].pos,
                              totals[l].pos + c.isize( ) );    }    }    }
          yreports[q] = xout.str( );    }
     int yc = 1;
     for ( int i = 0; i < yreports.isize( ); i++ )
     {    if ( yreports[i].size( ) > 0 ) 
          {    cout << "\n[" << yc++ << "] " << yreports[i];    }    }

     // Figure out the extent to which we've closed gaps.

     String gap_closure_summary;
     if (VALIDATE)
     {    vec< triple<int,int,int> > reftigs;
          int old_gaps = 0;
          for ( int i = 0; i < ptl_loc.isize( ); i++ )
          {    int j;
               for ( j = i + 1; j < ptl_loc.isize( ); j++ )
               {    if ( ptl_loc[j].first != ptl_loc[i].first ) break;
                    if ( ptl_loc[j-1].third - ptl_loc[j].second != K - 1 ) break;   }
               if ( reftigs.nonempty( ) 
                    && reftigs.back( ).first == ptl_loc[i].first )
               {    old_gaps++;    }
               reftigs.push( 
                    ptl_loc[i].first, ptl_loc[i].second, ptl_loc[j-1].third );
               i = j - 1;    }
          vec<Bool> join_to_next( reftigs.size( ), False );
          for ( int i = 0; i < newtigs.isize( ); i++ )
          {    const triple<int,int,int>& n = newtigs[i];
               int low = lower_bound( reftigs.begin( ), reftigs.end( ), n )
                    - reftigs.begin( );
               int high = upper_bound( reftigs.begin( ), reftigs.end( ), n )
                    - reftigs.begin( );
               while ( low > 0 && n.first == reftigs[low-1].first 
                    && n.second < reftigs[low-1].third )
               {    low--;    }
               while ( high < reftigs.isize( ) && n.first == reftigs[high].first
                    && n.third > reftigs[high].second )
               {    high++;    }
               for ( int j = low; j < high - 1; j++ )
               {    const triple<int,int,int>& t1 = reftigs[j];
                    const triple<int,int,int>& t2 = reftigs[j+1];
                    if ( t1.third - n.second >= K-1 && n.third - t2.second >= K-1 )
                         join_to_next[j] = True;    }    }
          cout << "\n";
          for ( int i = 0; i < reftigs.isize( ) - 1; i++ )
          {    if ( join_to_next[i] )
               {    cout << "see join from " << reftigs[i].first << "."
                         << reftigs[i].third << "-" << reftigs[i+1].second
                         << " (" << reftigs[i+1].second - reftigs[i].third + K - 1
                         << " kmers were missing)\n";    }    }
          gap_closure_summary = ToString(Sum(join_to_next)) + " of " 
               + ToString(old_gaps) + " gaps closed\n";    }

     // Save results.

     cout << Date( ) << ": join complete" << endl;

     // Rebuild and write (fake) unipaths.

     if (WRITE)  {    
       // Easy way - not true unipaths
       
       all_paths.Sort( );
       U.clear( );
       U.ReadAll(unibasefile);
       U.Append(all_paths);
       U.WriteAll( run_dir + "/" + OUT_HEAD + ".unibases.k"+ KS );   
     }


     // Print summary stats.

     cout << "\nSUMMARY\n" << "reads processed = " << reads_processed;
     if (VALIDATE)
     {    cout << ", truth found = "
               << simple_good + found_alls << " ("
               << PERCENT_RATIO( 3, simple_good + found_alls, reads_processed )
               << "), multiple solutions = " << found_alls << " ("
               << PERCENT_RATIO( 3, found_alls, reads_processed ) << ")\n";    }
     cout << gap_closure_summary;
     cout << "\ntime used = " << TimeSince(clock) << endl;
     _exit(0);    }
