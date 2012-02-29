///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixLocal.
//
// Needs AlignReads to have been run, and also SamplePairedReadDistributions.
//
// This is a very preliminary version of a program to edit an efasta assembly
// using read alignments to it.
//
// Scan assembly to locate weak spots, then locally reassemble.  The following
// steps are followed:
// 1. Look for positions providing a signal that something might be wrong.
// 2. Define windows around signal bases.
// 3. Identify reads that touch the windows.
// 4. Assemble the reads, yielding a graph.
// 5. Clean the graph.
// 6. Find all paths through the graph.
// 7. Score each path.
// 8. Pick winning paths.
// 9. Form into efasta.
//
// Known problems.
//
// 1. Presence of spurious edges.  In cases of tandem repeats (and perhaps in other
// cases), there is a tendency for the assembly graph to contain edges that are
// present only because of uncorrected sequencing errors.  There are various ways
// in which these errors might be corrected:
// (a) Running PreCorrect might help.
// (b) For each kmer x, we could test each position p on it by forming the four
// mutated kmers x_p:b, b in {A,C,G,T}.  For each b we would find all its instances
// in the corrected reads and form the list q_p:b of the associated quality scores
// at p.  Based on these data, under appropriate circumstances we might decide that
// some instances x_p:b are incorrect, and delete them from the graph.
// (c) FindErrors might be modified to do some version of (b).
//
// 2. Explosion in mapping reads.  In cases where the assembly graph contains
// enough alternatives (for example two loops may be enough), and a read contains
// many very low quality bases, mapping of the read to the graph may explode,
// because there are may be many equally likely paths for the read through the 
// graph.
//
// 3. Difficulty in scoring mappings to tandem repeats.  When a read is mapped to a
// tandem repeat, a fundamental question is the scoring of mappings that either
// terminate in a repeat copy, or leave the repeat.  Because the ends of reads tend
// to have low copy, this can be very difficult.
//
// 4. Overall computational performance.  This code has only been tested on 
// bacterial genomes.  On larger genomes the computational requirements would be
// much larger, both because the genomes are larger and because the rate of signal
// occurrence would generally be much higher.
//
// 5. Not checking for signal from pairing deviation.  This would be easy.
//
// 6. Need to sort alternatives by decreasing likelihood.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "FastIfstream.h"
#include "MainTools.h"
#include "Qualvector.h"
#include "feudal/QualNibbleVec.h"
#include "Superb.h"
#include "efasta/EfastaTools.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/AssemblyCleanupTools.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPath.h"
#include "paths/ReadLoc.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/Unipath.h"
#include "paths/FindErrorsCore.h"

void Coverage( const HyperKmerPath& h, const vec<tagged_rpint>& pathsdb,
     vec<double>& cov )
{    cov.resize( h.EdgeObjectCount( ) );
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    const KmerPath& p = h.EdgeObject(e);
          longlong c = 0;
          for ( int j = 0; j < p.NSegments( ); j++ )
          {    for ( longlong x = p.Start(j); x <= p.Stop(j); x++ )
               {    vec<longlong> con;
                    Contains( pathsdb, x, con );
                    c += con.size( );    }    }
          cov[e] = double(c) / double( p.KmerCount( ) );    }    }

void AddToPileup( const read_loc& rl, const basevector& b, const qualvector& q,
     const basevector& tig, vec<dumbcall>& calls )
{    align a;
     rl.GetAlign( a, b, tig );
     int p1 = a.pos1( ), p2 = a.pos2( );
     for ( int j = 0; j < a.Nblocks( ); j++ )
     {    if ( a.Gaps(j) > 0 )
          {    for ( int u = 0; u < a.Gaps(j); u++ )
                    calls[p2+u].base[4] += q[p1];
               p2 += a.Gaps(j);    }
          if ( a.Gaps(j) < 0 )
          {    for ( int u = 0; u < -a.Gaps(j); u++ )
                    calls[p2].base[5] += q[p1];
               p1 -= a.Gaps(j);    }
          for ( int x = 0; x < a.Lengths(j); x++ )
          {    calls[p2].base[ b[p1] ] += q[p1];
               ++p1; ++p2;    }    }    }

// AllSubpathTargets.  Find all paths E1-...-Em in a given HyperBasevector hb 
// that would be targets for gap-free alignment of a read of length N.  The 
// requirements for this are that
//     N <= len(E1-...-Em) and
//     N > len(E2-...-Em-1) + 1.
// Each returned path is given as a list of edge ids.

Bool AllSubpathTargets( const HyperBasevector& hb, const int N,
     vec< vec<int> >& subpaths )
{
     subpaths.clear( );
     vec<int> to_left, to_right;
     hb.ToLeft(to_left), hb.ToRight(to_right);
     int K = hb.K( );
     vec< vec<int> > partials;
     for ( int e = 0; e < hb.EdgeObjectCount( ); e++ )
     {    vec<int> p;
          p.push_back(e);
          partials.push_back(p);    }
     while( partials.nonempty( ) )
     {    vec<int> p = partials.back( );
          partials.pop_back( );
          int len1 = hb.EdgeLengthBases( p[0] );
          for ( int j = 1; j < p.isize( ); j++ )
               len1 += hb.EdgeLengthBases( p[j] ) - (K-1);
          Bool constraint1 = ( N <= len1 );
          Bool constraint2 = True;
          if ( p.size( ) >= 2 )
          {    int len2 = hb.EdgeLengthBases( p[1] );
               for ( int j = 2; j < p.isize( ) - 1; j++ )
                    len2 += hb.EdgeLengthBases( p[j] ) - (K-1);
               if ( !( N > len2 + 1 ) ) constraint2 = False;    }
          if ( !constraint2 ) continue;
          if (constraint1) subpaths.push_back(p);
          if ( subpaths.size( ) > 1000 ) return False;
          int w = to_right[ p.back( ) ];
          for ( int j = 0; j < hb.From(w).isize( ); j++ )
          {    vec<int> q = p;
               q.push_back( hb.EdgeObjectIndexByIndexFrom( w, j ) );
               partials.push_back(q);    }    }
     return True;    }

// Delete initial edges that do not contain the left flank, and terminal edges 
// that do no contain the right flank.

void DeleteIllegalTerminators( HyperKmerPath& h, const KmerBaseBroker& kbb,
     const basevector& left, const basevector& right )
{    while(1)
     {    vec<int> to_delete;
          for ( int v = 0; v < h.N( ); v++ )
          {    if ( h.Source(v) )
               {    for ( int j = 0; j < h.From(v).isize( ); j++ )
                    {    basevector b = kbb.Seq( h.EdgeObjectByIndexFrom( v, j ) );
                         if ( !b.ToString( ).Contains( left.ToString( ) ) )
                         {    to_delete.push_back( h.EdgeObjectIndexByIndexFrom( 
                                   v, j ) );    }    }    }
               if ( h.Sink(v) )
               {    for ( int j = 0; j < h.To(v).isize( ); j++ )
                    {    basevector b = kbb.Seq( h.EdgeObjectByIndexTo( v, j ) );
                         if ( !b.ToString( ).Contains( right.ToString( ) ) )
                         {    to_delete.push_back( h.EdgeObjectIndexByIndexTo( 
                                   v, j ) );    }    }    }    }
          UniqueSort(to_delete);
          h.DeleteEdges(to_delete);
          h.RemoveDeadEdgeObjects( );
          h.RemoveEdgelessVertices( );
          h.RemoveUnneededVertices( );
          if ( to_delete.empty( ) ) break;    }    }

void PickPath( 
     // inputs:
     const int low, const int high, const efasta& etiglet, const basevector& tig, 
     const int estart, const int estop, const HyperKmerPath& h, 
     const KmerBaseBroker& kbb, const vecbasevector& basesy, 
     const vecqualvector& qualsy, const int maxread, const vec< vec<int> >& paths,
     // heuristics:
     const int max_expand_to, const int q_junk_max, const int max_non_losers,
     const int max_paths,
     // logging:
     const int TIG, const Bool QLT, const String& tmp_dir, 
     const Bool SHOW_VOTES, const String& data_dir, ostream& rout,
     // output:
     vec< triple<int,int,String> >& replacements )
{
     int K = h.K( );
     vec<basevector> bpaths;
     rout << "\nThere are " << paths.size( ) << " paths:\n";
     for ( int m = 0; m < paths.isize( ); m++ )
     {    rout << "[" << m+1 << "] ";
          basevector b = kbb.Seq( h.EdgeObject( paths[m][0] ) );
          rout << "[";
          for ( int j = 0; j < paths[m].isize( ); j++ )
          {    if ( j > 0 ) rout << "-";
               rout << BaseAlpha( paths[m][j] );    
               if ( j > 0 )
               {    b.resize( b.isize( ) - (K-1) );
                    b = Cat( b, kbb.Seq( h.EdgeObject( paths[m][j] ) ) );    }    }
          bpaths.push_back(b);
          rout << "]\n";    }    
     
     // Now add in the assembly sequences.

     vec<basevector> A;
     if ( !etiglet.ExpandTo( A, max_expand_to ) ) rout << "expansion exploded\n";
     Bool first = True;
     int extras = 0;
     for ( int j = 0; j < A.isize( ); j++ )
     {    if ( Member( bpaths, A[j] ) ) continue;
          if (first) rout << "extra paths from assembly:\n";
          first = False;
          extras++;
          rout << "[" << bpaths.size( ) + 1 << "] " << A[j] << "\n";
          bpaths.push_back( A[j] );    }
     rout << "Found " << extras << " extra paths.\n";
     if ( bpaths.empty( ) ) { rout << "no path\n"; return; }
     if ( bpaths.solo( ) ) { rout << "unique path\n"; return; }
     if ( bpaths.isize( ) > max_paths )
     {    rout << "too many paths\n";
          return;    }

     // Assess sequences versus reference.

     if (QLT)
     {    Mkdir777(tmp_dir);
          String qlt_fasta = tmp_dir + "/x.QLT." + "fasta";
          {    Ofstream( out, qlt_fasta );
               basevector x;
               x.Print(out); // dummy
               for ( int j = 0; j < bpaths.isize( ); j++ )
               {    const int big_flank = 8000;
                    int nleft = Min( big_flank, low );
                    int nright = Min( big_flank, tig.isize( ) - high );
                    basevector left( tig, low - nleft, nleft );
                    basevector right( tig, high, nright );
                    basevector b = Cat( left, bpaths[j], right );
                    b.Print( out, j+1 );    }    }
          SystemSucceed( "QueryLookupTable K=12 MM=12 MC=0.15 SEQS=" + qlt_fasta 
               + " L=" + data_dir + "/genome.lookup VISUAL=True PARSEABLE=True "
               "QUIET=True NH=True > " + tmp_dir + "/x.QLT.out" );
          vec<look_align> aligns;
          LoadLookAligns( tmp_dir + "/x.QLT.out", aligns );
          vec<int> errs( aligns.size( ) );
          for ( int j = 0; j < aligns.isize( ); j++ )
               errs[j] = aligns[j].Errors( );
          vec<int> ids( aligns.size( ), vec<int>::IDENTITY );
          SortSync( errs, ids );
          Bool winner_known = False;
          if ( errs.size( ) == 0 
               || ( errs.size( ) >= 2 && errs[0] == errs[1] )
               || !aligns[ ids[0] ].FullLength( ) )
          {    rout << "\nWarning: no clear winner\n";    }
          else
          {    winner_known = True;
               int w = ids[0];
               rout << "\nReference says winner is [" << w+1 
                    << "], errors = " << errs[0] << "\n";    }
          if ( !winner_known )
          {    rout << "\nAlignments:\n";
               fast_ifstream in( tmp_dir + "/x.QLT.out" );
               String line;
               while(1)
               {    getline( in, line );
                    if ( in.fail( ) ) break;
                    if ( !line.Contains( "QUERY", 0 ) ) 
                         rout << line << "\n";    }    }
          SystemSucceed( "/bin/rm -rf " + tmp_dir );    }

     // Extend paths using flanking sequence.

     vec<basevector> bpaths_orig(bpaths);
     int left_flank_size = Min( maxread, low );
     int right_flank_size = Min( maxread, tig.isize( ) - high );
     basevector left_flank( tig, low - left_flank_size, left_flank_size );
     basevector right_flank( tig, high, right_flank_size );
     for ( int j = 0; j < bpaths.isize( ); j++ )
          bpaths[j] = Cat( left_flank, bpaths[j], right_flank );

     // Align the reads to the extended paths, or "targets".  Note 
     // the following issues:
     // (1) We do not allow indels.
     // (2) We do not take into account consistency with the partner placement.
     //     Thus the best placement for the read taken individually might not be the
     //     best placement overall.
     // Voting.  For each gap-free placement of a read, we form the sum of its
     // quality scores at mismatches.  Quality scores of two are less are treated as
     // zero.  We find the target having the lowest score (the "winner"), and let
     // delta be the difference between the score of the runner up and the winner.
     // The target is then assigned a vote of 1 - 10^(-delta/10).

     vec< vec<double> > VOTE( bpaths.size( ) );
     for ( int j = 0; j < bpaths.isize( ); j++ )
          VOTE[j].resize( bpaths.size( ), 0 );
     vec< vec<int> > SCORE( basesy.size( ) );
     for ( size_t j = 0; j < basesy.size( ); j++ )
     {    SCORE[j].resize( bpaths.size( ), 1000000000 );
          for ( int i = 0; i < bpaths.isize( ); i++ )
          {    const basevector &bp = bpaths[i], &r = basesy[j];
               const qualvector& q = qualsy[j];
               for ( int start = 0; start <= bp.isize( ) - r.isize( ); start++ )
               {    int qsum = 0;
                    for ( int l = 0; l < r.isize( ); l++ )
                    {    if ( r[l] != bp[start+l] )
                         {    double qx = q[l];
                              if ( qx <= q_junk_max ) qx = 0;
                              qsum += qx;
                              if ( qsum >= SCORE[j][i] ) break;    }    }
                    SCORE[j][i] = Min( SCORE[j][i], qsum );    }    }    }
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     for ( int i2 = i1 + 1; i2 < bpaths.isize( ); i2++ )
     {    vec<basevector> BPATHS;
          BPATHS.push_back( bpaths[i1], bpaths[i2] );
          vec<double> vote( 2, 0 );
          for ( size_t j = 0; j < basesy.size( ); j++ )
          {    vec<int> score(2);
               score[0] = SCORE[j][i1], score[1] = SCORE[j][i2];
               vec<int> ids( BPATHS.size( ), vec<int>::IDENTITY );
               /*
               rout << "scores for read " << j; // XXXXXXXXXXXXXXXXXXXXXXX
               for ( int m = 0; m < score.isize( ); m++ ) // XXXXXXXXXXXXX
                    rout << " " << score[m]; // XXXXXXXXXXXXXXXXXXXXXXXXXX
               rout << "\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               */
               SortSync( score, ids );
               int delta = score[1] - score[0];
               if ( delta > 0 )
                    vote[ ids[0] ] += 1.0 - pow( 10.0, -double(delta)/10.0 );    }
          VOTE[i1][i2] = vote[0], VOTE[i2][i1] = vote[1];
          if (SHOW_VOTES)
          {    rout << "votes:\n";
               for ( int j = 0; j < BPATHS.isize( ); j++ )
               {    int id = ( j == 0 ? i1+1 : i2+1 );
                    rout << "[" << id << "] " << vote[j] 
                         << "\n";    }    }    }
     
     int winner = -1;
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    Bool i1_wins = True;
          for ( int i2 = 0; i2 < bpaths.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               if ( VOTE[i1][i2] == 0 || VOTE[i2][i1] >= 2.0 
                    || VOTE[i1][i2] < 10.0 * VOTE[i2][i1] )
               {    i1_wins = False;    }    }
          if (i1_wins) winner = i1;    }
     if ( winner >= 0 )
     {    rout << "winner: [" << winner+1 << "]\n";
          if ( bpaths_orig[winner].ToString( ) == etiglet )
               rout << "equals assembly\n";
          else
          {    rout << "replacing assembly sequence\n" << etiglet 
                    << "\nby\n" << bpaths_orig[winner].ToString( ) << "\n";
               PRINT3_TO( rout, TIG, estart, estop );
               replacements.push( estart, estop,
                    bpaths_orig[winner].ToString( ) );    }
          return;    }

     vec<Bool> loser( bpaths.size( ), False );
     vec<Bool> beats_something( bpaths.size( ), False );
     vec<Bool> something_beats( bpaths.size( ), False );

     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    for ( int i2 = 0; i2 < bpaths.isize( ); i2++ )
          {    if ( i1 == i2 ) continue;
               if ( VOTE[i2][i1] < 2.0 && VOTE[i1][i2] >= 10.0 * VOTE[i2][i1] )
               {    beats_something[i1] = True;
                    something_beats[i2] = True;    }    }    }
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
          if ( something_beats[i1] && !beats_something[i1] ) loser[i1] = True;
     rout << "non-losers:\n";
     int n_non_losers = 0;
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     {    if ( !loser[i1] ) 
          {    rout << "[" << i1+1 << "]\n";
               n_non_losers++;    }    }
     vec<int> vote_sum( bpaths.size( ), 0 );
     for ( int i1 = 0; i1 < bpaths.isize( ); i1++ )
     for ( int i2 = i1 + 1; i2 < bpaths.isize( ); i2++ )
     {    if ( !loser[i1] && !loser[i2] )
          {    rout << "vote:\n" << "[" << i1+1 << "] " << VOTE[i1][i2] << "\n";
               rout << "[" << i2+1 << "] " << VOTE[i2][i1] << "\n";    
               vote_sum[i1] += VOTE[i1][i2], vote_sum[i2] += VOTE[i2][i1];    }    }
     if ( n_non_losers <= max_non_losers )
     {    vec<basevector> nl;
          vec<int> nl_votes;
          for ( int j = 0; j < bpaths.isize( ); j++ )
          {    if ( !loser[j] ) 
               {    nl.push_back( bpaths_orig[j] );
                    nl_votes.push_back( vote_sum[j] );    }    }
          ReverseSortSync( nl_votes, nl );
          if ( nl != A )
          {    rout << "replacing assembly sequence\n" << etiglet 
                    << "\nby\n" << efasta(nl) << "\n";
               PRINT3_TO( rout, TIG, estart, estop );
               replacements.push( estart, estop, efasta(nl) );    }    }    }

class heuristics {

     public:

     int64_t max_reads;
     int min_calls;
     double agree_ceil_weak;
     double agree_floor_strong;
     double MIN_FRAC;
     int q_junk_max;
     int max_tandem_period;
     int min_tandem_copies;
     int min_tandem_length;
     double max_offby;
     int max_non_losers;
     int max_expand_to;
     int max_paths;
     double gapdev_mult;
     int boundary_push;
     int max_cyclic_paths;
     int max_cyclic_loops;
     double min_devs_off;
     Bool do_cyclic;
     int max_cyclic_edges;

};

class logging_control {

     public:

     String DUMP_EC;
     Bool SHOW_VOTES;
     Bool QLT;
     String data_dir;
     String DOT;

};

const int L = 40;

void ProcessWindow(
     // inputs:
     const vec< triple<kmer<L>,int,int> >& kmers_plus, const vec< kmer<L> >& kmers,
     const heuristics& heur, const logging_control& logc, const String& tmp_dir,
     pair<int,int>& window, const String& signal, const basevector& tig,
     const efasta& T, const int TIG, const int K, const vec<int>& ids_frag,
     const vec<int>& ids_jump, const vec<read_loc>& locs,
     const vec<int>& loc_id_frag, const vec<int>& loc_id_jump,
     const vecbasevector& bases, const vecqualvector& quals,
     const vec<int>& readlengths, const vecbasevector& tigs,
     const vec< vec<longlong> >& aligns_index, const vec<Bool>& placed_fw, 
     const vec<Bool>& placed_rc, const vec< pair<int,int> >& placement,
     const int nlibs, const int nlibs_frag, const vec< vec<int> >& DISTS,
     const vec< vec<double> >& X, const vec<Bool>& lfail,
     const vec<unsigned short>& read_lengths, const PairsManager& pairs,
     // outputs:
     ostream& rout, vec< triple<int,int,int> >& markups,
     vec< triple<int,int,String> >& replacements )
{
     rout << signal << "\n";
     restart:
     int low = window.first, high = window.second; 

     // Set up for PredictGap computation.  For now, to avoid unneeded calls, we 
     // don't do it here.

     int GAP;
     double DEV;
     const int GVERBOSITY = 0;
     const int max_overlap = 0;
     const vec<int> libs_to_use;

     // Get the contig chunk and its flanks.

     rout << "\n";
     basevector left( tig, low, K ), right( tig, high - K, K );
     basevector tiglet( tig, low, high - low );
     int pos = 0, estart = -1, estop = -1;
     for ( int j = 0; j < T.isize( ); j++ )
     {    if ( pos == low ) estart = j;
          if ( pos == high ) estop = j;
          if ( T[j] == '{' )
          {    for ( j++; j < T.isize( ); j++ )
               {    if ( T[j] == ',' ) break;
                    pos++;    }
               for ( j++; j < T.isize( ); j++ )
                    if ( T[j] == '}' ) break;    }
          else pos++;    }
     efasta etiglet = efasta( T.substr( estart, estop - estart ) );
     etiglet.Print( rout, "etiglet (" + ToString( etiglet.Length1( ) ) + " bases)" );

     // Isolate the reads for this region.

     vecbasevector basesy;
     vecqualvector qualsy;
     for ( int j = 0; j < ids_frag.isize( ); j++ )
     {    const read_loc& rl = locs[ loc_id_frag[j] ];
          if ( IntervalOverlap( rl.Start( ), rl.Stop( ), low, high ) > 0 )
          {    basesy.push_back( bases[j] );
               qualsy.push_back( quals[j] );    }    }
     for ( int j = 0; j < ids_jump.isize( ); j++ )
     {    const read_loc& rl = locs[ loc_id_jump[j] ];
          if ( IntervalOverlap( rl.Start( ), rl.Stop( ), low, high ) > 0 )
          {    basesy.push_back( bases[ j + ids_frag.isize( ) ] );
               qualsy.push_back( 
                    quals[ j + ids_frag.isize( ) ] );    }    }
     rout << "have " << basesy.size( ) << " reads" << endl;
     if ( (int64_t) basesy.size( ) > heur.max_reads )
     {    rout << "too many reads, giving up" << endl;
          return;    }

     // Correct errors in the reads.

     BaseVecVec basesx = basesy;
     {    const size_t nqv = qualsy.size();
          QualNibbleVecVec qualsx(nqv);
          for (size_t iqv = 0; iqv < nqv; iqv++)
          {    for (size_t iq = 0; iq < qualsy[iqv].size(); iq++)
               {    qualsx[iqv].push_back(qualsy[iqv][iq]);    }    }
          const unsigned K = 24;
          const unsigned n_cycles = 2;
          const unsigned verbosity = 0;
          find_errors(efp_default, K, &basesx, &qualsx, n_cycles, verbosity);
          if ( logc.DUMP_EC != "" ) basesx.WriteAll(logc.DUMP_EC);    }

     vecKmerPath Paths, Paths_rc, unipaths;
     vec<tagged_rpint> Pathsdb, unipathsdb;
     ReadsToPathsCoreY( basesx, K, Paths );
     CreateDatabase( Paths, Paths_rc, Pathsdb );
     Unipath( Paths, Paths_rc, Pathsdb, unipaths, unipathsdb );
     KmerBaseBroker kbb( K, Paths, Paths_rc, Pathsdb, basesx );
     digraph AA;
     BuildUnipathAdjacencyGraph(Paths, Paths_rc, Pathsdb, unipaths, unipathsdb, AA);
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, AA, unipaths, h );
     int MIN_COMPONENT = 20;
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     rout << "initial graph has " << h.EdgeObjectCount( ) << " edges" << endl;

     // Delete initial edges that do not contain the left flank, and
     // terminal edges that do no contain the right flank.

     DeleteIllegalTerminators( h, kbb, left, right );

     // If an initial [resp. terminal] edge contains the left [resp. right]
     // flank, truncate it there.

     for ( int v = 0; v < h.N( ); v++ )
     {    if ( h.Source(v) )
          {    for ( int j = 0; j < h.From(v).isize( ); j++ )
               {    KmerPath& e = h.EdgeObjectByIndexFromMutable( v, j );
                    basevector b = kbb.Seq(e);
                    int p = b.ToString( ).Position( left.ToString( ) );
                    if ( p > 0 )
                    {    KmerPath enew;
                         e.CopySubpath( e.Begin( ) + p, e.End( ), enew );
                         e = enew;    }    }    }
          if ( h.Sink(v) )
          {    for ( int j = 0; j < h.To(v).isize( ); j++ )
               {    KmerPath& e = h.EdgeObjectByIndexToMutable( v, j );
                    basevector b = kbb.Seq(e);
                    int p = b.ToString( ).Position( right.ToString( ) );
                    if ( p >= 0 && p + K < b.isize( ) )
                    {    KmerPath enew;
                         e.CopySubpath( e.Begin( ), 
                              e.End( ) - ( b.isize( ) - p - K), enew );
                         e = enew;    }    }    }    }

     // Prune branches using MIN_FRAC.

     vec<double> cov;
     Coverage( h, Pathsdb, cov );
     vec<int> to_delete;
     for ( int v = 0; v < h.N( ); v++ )
     {    for ( int j1 = 0; j1 < h.From(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.From(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexFrom( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexFrom( v, j2 );
                    // if ( h.From(v)[j1] != h.From(v)[j2] ) continue;
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 >= 2.0 && c2 >= 2.0 ) continue;
                    if ( c1 < int( floor( heur.MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }
          for ( int j1 = 0; j1 < h.To(v).isize( ); j1++ )
          {    for ( int j2 = 0; j2 < h.To(v).isize( ); j2++ )
               {    int e1 = h.EdgeObjectIndexByIndexTo( v, j1 );
                    int e2 = h.EdgeObjectIndexByIndexTo( v, j2 );
                    // if ( h.To(v)[j1] != h.To(v)[j2] ) continue;
                    double c1 = cov[e1], c2 = cov[e2];
                    if ( c1 >= 2.0 && c2 >= 2.0 ) continue;
                    if ( c1 < int( floor( heur.MIN_FRAC * (c1+c2) ) ) )
                         to_delete.push_back(e1);    }    }    }
     UniqueSort(to_delete);
     h.DeleteEdges(to_delete);
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Look for grubby kmers.

     vec<int> grubby;
     const int KK(20);
     vec< kmer<KK> > B, BX;
     vec<qualvector> Q;
     kmer<KK> x;
     int xcount = 0;
     for ( size_t j = 0; j < basesy.size( ); j++ )
          xcount += basesy[j].isize( ) - KK + 1;
     B.reserve(xcount), Q.reserve(xcount);
     for ( size_t j = 0; j < basesy.size( ); j++ )
     {    for ( int l = 0; l <= basesy[j].isize( ) - KK; l++ )
          {    x.SetToSubOf( basesy[j], l );
               qualvector q;
               q.SetToSubOf( qualsy[j], l, KK );
               B.push_back(x), Q.push_back(q);    }    }
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    const basevector& b = kbb.Seq( h.EdgeObject(j) );
          for ( int l = 0; l <= b.isize( ) - KK; l++ )
          {    x.SetToSubOf( b, l );
               BX.push_back(x);    }    }
     SortSync( B, Q ), Sort(BX);
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    basevector b = kbb.Seq( h.EdgeObject(j) );
          for ( int l = 0; l <= b.isize( ) - KK; l++ )
          {    basevector c( b, l, KK );
               for ( int r = 0; r < KK; r++ )
               {    int m0 = c[r];
                    vec<Bool> ref( 4, False );
                    for ( int m = 0; m < 4; m++ )
                    {    if ( m == m0 ) ref[m] = True;
                         else 
                         {    basevector d(c);
                              d.Set( r, m );
                              x.Set(d);
                              if ( BinMember( BX, x ) ) ref[m] = True;    }    }
                    if ( Sum(ref) == 1 ) continue;
                    vec<int> low(4), high(4);
                    for ( int m = 0; m < 4; m++ )
                    {    basevector d(c);
                         d.Set( r, m );
                         x.Set(d);
                         low[m] = LowerBound( B, x );
                         high[m] = UpperBound( B, x );    }
                    Bool have_nonref = False;
                    for ( int m = 0; m < 4; m++ )
                    {    if ( m == m0 ) continue;
                         if ( low[m] < high[m] ) have_nonref = True;    }
                    if ( !have_nonref ) continue;

                    // Kill the kmer if it is a mutation of another kmer,
                    // if the other kmer has at least 10 times as much 
                    // support, and if the maximum quality score at the 
                    // mutation position is < the median quality score of 
                    // the other kmer.
     
                    const int cov_qtest_ratio = 10;
                    vec< vec<int> > qs(4);
                    for ( int m = 0; m < 4; m++ )
                    {    for ( int x = low[m]; x < high[m]; x++ )
                              qs[m].push_back( Q[x][r] );
                         Sort( qs[m] );    }
                    for ( int m = 0; m < 4; m++ )
                    {    if ( qs[m0].size( ) * cov_qtest_ratio 
                              > qs[m].size( ) )
                         {    continue;    }
                         if ( qs[m0].nonempty( ) && qs[m0].back( ) <
                              Median( qs[m] ) )
                         {    grubby.push_back(j);    }    }
                    Bool ec_verbose = False;
                    if (ec_verbose)
                    {    rout << "\n" << c << " " << r << "\n";
                         for ( int m = 0; m < 4; m++ )
                         {    const vec<int>& v = qs[m];
                              rout << as_base(m) << (m == c[r] ? " (ref" : "");
                              if ( v.isize( ) < 20 )
                              {    for ( int x = 0; x < v.isize( ); x++ )
                                        rout << " " << v[x];    }
                              else
                              {    rout << " " << v.size( ) << " quals of median " 
                                        << v[ v.size( )/2 ];    }
                              rout << "\n";    }    }    }    }    }
     UniqueSort(grubby);
     h.DeleteEdges(grubby);
     if ( MIN_COMPONENT > 0 ) h.RemoveSmallComponents(MIN_COMPONENT);
     h.RemoveDeadEdgeObjects( );
     h.RemoveEdgelessVertices( );
     h.RemoveUnneededVertices( );
     h.RemoveDeadEdgeObjects( );

     // Delete initial edges that do not contain the left flank, and
     // terminal edges that do no contain the right flank.

     DeleteIllegalTerminators( h, kbb, left, right );

     // If left and right are not both present, kill the graph.

     int left_id = -1, right_id = -1;
     vec<int> to_left, to_right;
     h.ToLeft(to_left), h.ToRight(to_right);
     for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
     {    if ( kbb.Seq( h.EdgeObject(j) ).ToString( ).Contains( left.ToString( ) ) )
          {    left_id = to_left[j];    }
          if ( kbb.Seq( h.EdgeObject(j) ).ToString( ).Contains( right.ToString( ) ) )
          {    right_id = to_right[j];    }    }
     if ( left_id < 0 || right_id < 0 )
     {    rout << "graph is empty" << endl;
          if ( !etiglet.Contains( "{" ) )
          {    rout << "marking assembly sequence\n";
               PRINT3_TO( rout, TIG, estart, estop );
               markups.push( TIG, estart, estop );    }
          return;    }

     // Present results.

     rout << "\n";
     PRINT_TO( rout, h.EdgeObjectCount( ) );
     if ( logc.DOT != "" )
     {    Ofstream( out, logc.DOT );
          h.PrintSummaryDOT0w( out, False, False, True, NULL, False );    }
     Coverage( h, Pathsdb, cov );
     for ( int e = 0; e < h.EdgeObjectCount( ); e++ )
     {    int v = to_left[e], w = to_right[e];
          rout << ">edge_" << BaseAlpha(e) << " " << v << ":" << w 
               << " kmers=" << h.EdgeLengthKmers(e)
               << " cov=" << ToString( cov[e] ) << "\n";
          kbb.ToSequence( h.EdgeObject(e) ).PrintN(rout);    }    

     // Handle the case where there are cycles.

     if ( !h.Acyclic( ) )
     {    rout << "cyclic\n";

          // Give up if too many edges.

          if ( h.EdgeObjectCount( ) > heur.max_cyclic_edges )
          {    rout << "\ntoo many edges in cyclic graph, giving up\n";
               return;    }

          // To proceed we require that the graph has a unique source and
          // sink, and unique edges emanating from them.  If we don't
          // find a source or a sink, we extend the boundaries of the
          // window, which is done via a goto.  Yecch.

          vec<int> sources, sinks;
          h.Sources(sources), h.Sinks(sinks);
          PRINT2_TO( rout, sources.size( ), sinks.size( ) );
          if ( sources.empty( ) )
          {    if ( low == 0 )
               {    rout << "\nno sources, can't push left boundary "
                         << "further, giving up\n";
                    return;    }
               low -= heur.boundary_push;
               if ( low < 0 ) low = 0;

               // Lower low if it lands illegally.

               int e_low = T.Index1Alt(low);
               int brack = 0;
               for ( int j = 0; j < e_low; j++ )
               {    if ( T[j] == '{' ) brack++;
                    if ( T[j] == '}' ) brack--;    }
               if ( brack > 0 )
               {    while ( low > 0 )
                    {   low--;
                        e_low--;
                        if ( T[e_low] == '{' ) break;    }    }

               // Push boundary.

               window.first = low;
               rout << "\n*** no sources, pushing left boundary of "
                    << "window to " << low << " ***\n\n";
               goto restart;    }
          if ( sinks.empty( ) )
          {    if ( high == tig.isize( ) )
               {    rout << "\nno sinks, can't push right boundary "
                         << "further, giving up\n";
                    return;    }
               high += heur.boundary_push;
               if ( high > tig.isize( ) ) high = tig.size( );

               // Raise high if it lands illegally.

               get_high:
               int e_high = T.Index1Alt(high);
               int brack = 0;
               for ( int j = 0; j < e_high; j++ )
               {    if ( T[j] == '{' ) brack++;
                    if ( T[j] == '}' ) brack--;    }
               if ( brack > 0 && high < tig.isize( ) )
               {    high++;
                    goto get_high;    }

               // Push boundary.

               window.second = high;
               rout << "\n*** no sinks, pushing right boundary of window "
                    << "to " << high << " ***\n\n";
               goto restart;    }
          if ( !sources.solo( ) || !sinks.solo( ) )
          {    rout << "don't have unique source and sink\n";
               return;    }
          int v = sources[0], w = sinks[0];
          int source_edges = h.From(v).size( );
          int sink_edges = h.To(w).size( );
          if ( source_edges != 1 || sink_edges != 1 )
          {    rout << "have unique source and sink\n";
               PRINT2_TO( rout, source_edges, sink_edges );
               rout << "don't have unique edges emanating from "
                    << "source and sink\n";
               return;    }
          rout << "have unique source and sink and edges "
               << "emanating from them\n";
          rout << BaseAlpha( h.EdgeObjectIndexByIndexFrom( v, 0 ) )
               << " --> ... --> "
               << BaseAlpha( h.EdgeObjectIndexByIndexTo( w, 0 ) ) << "\n";
          HyperBasevector hb( h, kbb );
          int in_id = hb.EdgeObjectIndexByIndexFrom( v, 0 );
          int out_id = hb.EdgeObjectIndexByIndexTo( w, 0 );
          if ( !hb.EdgeObject(in_id).ToString( ).Contains( 
               left.ToString( ), 0 )
               || !hb.EdgeObject(out_id).ToString( ).Contains( 
               right.ToString( ), -1 ) )
          {    rout << "can't find left and right\n";
               return;    }

          // Compute gap.

          Bool GAP_PREDICTED;
          PredictGap( "jump", kmers_plus, kmers, max_overlap, nlibs, nlibs_frag, 
               libs_to_use, lfail, tigs, DISTS, X, bases, pairs, read_lengths, 
               aligns_index, placed_fw, placed_rc, placement, 
               gap_id( gap_id::WITHIN, TIG, window.first, window.second ), 
               GVERBOSITY, GAP_PREDICTED, GAP, DEV );
          double devs_off = 0.0;
          Bool way_off = False;
          if (GAP_PREDICTED)
          {    GAP = ( window.second - window.first ) - GAP;
               PRINT2_TO( rout, GAP, DEV );
               devs_off = Abs( double(GAP)/DEV );
               way_off = ( devs_off >= heur.min_devs_off );    }
          else rout << "Failed to predict gap.\n";

          // Could there be only a small number of paths through the
          // graph that are consistent with distance constraints?

          if (GAP_PREDICTED)
          {    vec<int> LL;
               for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
                    LL.push_back( h.EdgeLengthKmers(j) );
               digraphE<int> G( h, LL );
               int ex = etiglet.Length1( ) - (K-1);
               int LL1 = Max( 0, ex - int( ceil( double(GAP) 
                    + heur.gapdev_mult * DEV ) ) );
               int LL2 = Max( 0, ex - int( floor( double(GAP) 
                    - heur.gapdev_mult * DEV ) ) );
               vec< vec<int> > paths;
               if ( G.AllPathsLengthRange( v, w, LL1, LL2, to_right, paths,
                    heur.max_cyclic_paths, heur.max_cyclic_loops ) )
               {    rout << "\nsee " << paths.size( ) 
                         << " paths through the graph that are "
                         << "consistent with distance constraints\n";    
                    if ( paths.nonempty( ) )
                    {    PickPath( low, high, etiglet, tig, estart, estop, h, kbb, 
                              basesy, qualsy, Max(readlengths), paths, 
                              heur.max_expand_to, heur.q_junk_max, 
                              heur.max_non_losers, heur.max_paths,
                              TIG, logc.QLT, tmp_dir, logc.SHOW_VOTES, 
                              logc.data_dir, rout, replacements );    
                         return;    }    }    }

          // Test for simple cyclic graph.

          Bool simple_cyclic = False;
          if ( h.EdgeObjectCount( ) == 4 )
          {    if ( h.From(v).solo( ) && h.To(w).solo( ) )
               {    int x = h.From(v)[0], y = h.To(w)[0];
                    if ( h.From(x).solo( ) && h.To(x).size( ) == 2
                         && h.From(y).size( ) ==2 && h.To(y).solo( ) )
                    {    simple_cyclic = True;
                         rout << "simple cyclic\n";    }    }    }
          if ( !simple_cyclic && !heur.do_cyclic ) 
          {    if (way_off) 
               {    rout << "way off\n" << "marking assembly sequence\n";
                    PRINT3_TO( rout, TIG, estart, estop );
                    markups.push( TIG, estart, estop );    }
               rout << "giving up, as the graph is cyclic but not simple" << endl;
               return;    }
     
          // Create an extended version of the graph in which we add
          // a readlength of bases to both ends.

          HyperBasevector hbplus(hb);
          basevector& in_edge = hbplus.EdgeObjectMutable(in_id);
          basevector& out_edge = hbplus.EdgeObjectMutable(out_id);
          int mr = Max(readlengths);
          int lext = Min( mr, low ), rext = Min( mr, tig.isize( ) - high );
          basevector left_ext( tig, low - lext, lext );
          basevector right_ext( tig, high, rext );
          in_edge = Cat( left_ext, in_edge );
          out_edge = Cat( out_edge, right_ext );

          // Align the reads to the edges.

          vec<int> source_sink_counts;
          Bool exploded = False;
          vec< vec< vec<int> > > subpaths( readlengths.size( ) );
          for ( int r = 0; r < readlengths.isize( ); r++ )
          {    exploded = !AllSubpathTargets( 
                    hbplus, readlengths[r], subpaths[r] );
               if (exploded)
                    rout << "Warning, AllSubpathTargets exploded\n";    }
          for ( size_t r = 0; r < basesy.size( ); r++ )
          {    const basevector& R = basesy[r];
               const qualvector& Q = qualsy[r];
               int px = BinPosition( readlengths, R.isize( ) );
               vec<int> score( subpaths[px].size( ) );
               for ( int j = 0; j < subpaths[px].isize( ); j++ )
               {    const vec<int>& p = subpaths[px][j];
                    basevector target = hbplus.EdgeObject( p[0] );
                    for ( int l = 1; l < p.isize( ); l++ )
                    {    target.resize( target.isize( ) - (K-1) );
                         target = Cat( 
                              target, hbplus.EdgeObject( p[l] ) );    }
                    int n1 = hbplus.EdgeLengthBases( p[0] );
                    int bestq = 1000000000;
                    int nx = n1;
                    for ( int l = 1; l < p.isize( ) - 1; l++ )
                         nx += hbplus.EdgeLengthKmers( p[l] );
                    for ( int start = Max( 0, nx - R.isize( ) + 1 ); 
                         start < n1 - (K-1); start++ )
                    {    if ( start + R.isize( ) - 1 < target.isize( ) )
                         {    int qsum = 0;
                              for ( int l = 0; l < R.isize( ); l++ )
                              {    if ( R[l] != target[start+l] ) 
                                        qsum += Q[l];    }
                              bestq = Min( bestq, qsum );    }    }
                    score[j] = bestq;    }    
               vec<int> ids( subpaths[px].size( ), vec<int>::IDENTITY );
               SortSync( score, ids );
               const int max_diff = 20;
               const int max_score = 200;
               if ( score[0] <= max_score )
               {    rout << "\nalignments of read " << r << "\n";
                    rout << R.ToString( ) << "\n";
                    Bool source_sink = True;
                    vec< vec<String> > rows;
                    for ( int j = 0; j < score.isize( ); j++ )
                    {    if ( score[j] > score[0] + max_diff ) break;
                         vec<String> row;
                         row.push_back( "[" + ToString(ids[j]+1) + "]" );
                         row.push_back( "(score=" + ToString(score[j]) + ")" );   
                         const vec<int>& p = subpaths[px][ ids[j] ];
                         String r;
                         for ( int l = 0; l < p.isize( ); l++ )
                         {    r += ( l > 0 ? "-" : "" ) 
                                   + BaseAlpha( p[l] );    }
                         row.push_back(r);
                         rows.push_back(row);
                         if ( !h.Source( to_left[ p.front( ) ] )
                              || !h.Sink( to_right[ p.back( ) ] ) )
                         {    source_sink = False;    }    }
                    PrintTabular( rout, rows, 1, "lll" );
                    if (source_sink)
                    {    source_sink_counts.push_back(
                              ( subpaths[px][ ids[0] ].isize( ) - 3 ) / 2 );
                         rout << "This read aligns as a "
                              << "source-sink.\n";    }    }    }
          rout << "\n";
          PRINT_TO( rout, source_sink_counts.size( ) );
          Sort(source_sink_counts);
          for ( int i = 0; i < source_sink_counts.isize( ); i++ )
          {    int j = source_sink_counts.NextDiff(i);
               rout << "source-sink(" << source_sink_counts[i] << " copies): "
                    << j - i << " times\n";
               i = j - 1;    }
          if (GAP_PREDICTED) PRINT2_TO( rout, GAP, DEV );
                              
          // Unwind etiglet along graph.

          int x = h.From(v)[0];
          int fwd_id = h.EdgeObjectIndexByIndexFrom( x, 0 );
          int y = h.From(x)[0];
          int rev_id = -1;
          for ( int j = 0; j < h.From(y).isize( ); j++ )
          {    if ( h.From(y)[j] == x )
                    rev_id = h.EdgeObjectIndexByIndexFrom( y, j );    }
          int xcount = -1;
          for ( int count = 0; ; count++ )
          {    vec<int> v;
               v.push_back(in_id);
               for ( int r = 0; r < count; r++ )
                    v.push_back( fwd_id, rev_id );
               v.push_back( fwd_id, out_id );
               basevector e = hb.EdgePathToBases(v);
               if ( e.isize( ) > etiglet.isize( ) ) break;
               if ( e.ToString( ) == etiglet )
               {    xcount = count;
                    break;    }    }
          int repeat_size = h.EdgeLengthKmers(fwd_id) + h.EdgeLengthKmers(rev_id);
          if ( xcount >= 0 )
          {    rout << "observe " << xcount << " repeat copies in existing "
                    << "contig, each of size " << repeat_size << "\n";    }
          else 
          {    rout << "can't map sequence in contig to graph\n";
               if ( !etiglet.Contains( "{" ) )
               {    rout << "marking assembly sequence\n";
                    PRINT3_TO( rout, TIG, estart, estop );
                    markups.push( TIG, estart, estop );    }
               return;   }

          // Decide what to do next.

          const int min_source_sinks = 3;
          if ( source_sink_counts.isize( ) < min_source_sinks && GAP_PREDICTED )
          {    rout /* << "Can't find enough source-sink edges, suggest " */
                    << "using GAP = " << GAP << " +/- " 
                    << DEV << " to tweak\nnumber of copies of repeat.\n";
               if ( simple_cyclic && !exploded /* && nbounds > 0 */
                    && !signal.Contains( "ambiguity" ) )
               {    rout << "to perturb tandem count\n";    
                    int add_low = int( floor( ( double(-GAP) 
                         - heur.gapdev_mult * DEV ) / double(repeat_size) ) );
                    int add_high = int( ceil( ( double(-GAP) 
                         + heur.gapdev_mult * DEV ) / double(repeat_size) ) );
                    rout << "to add between " << add_low << " and "
                         << add_high << " copies of repeat\n";    
                    int new_count_low = Max( 0, xcount + add_low );
                    int new_count_high = Max( 0, xcount + add_high );
                    int expected = Max( 0, add_low 
                         + int( round( double(-GAP)/double(repeat_size) ) ) );
                    PRINT3_TO( rout, new_count_low, expected, new_count_high );
                    if ( !( new_count_low <= new_count_high ) )
                    {    rout << "doesn't make sense\n";
                         return;    }

                    // Replace tandem repeat.

                    vec<basevector> answer;
                    vec<int> counts;
                    counts.push_back(expected);
                    for ( int c = new_count_low; c <= new_count_high; c++ )
                         if ( c != expected ) counts.push_back(c);
                    for ( int ci = 0; ci < counts.isize( ); ci++ )
                    {    int c = counts[ci];
                         vec<int> v;
                         v.push_back(in_id);
                         for ( int r = 0; r < c; r++ )
                              v.push_back( fwd_id, rev_id );
                         v.push_back( fwd_id, out_id );
                         answer.push_back( hb.EdgePathToBases(v) );    }
                    rout << "replacing assembly sequence\n" 
                         << etiglet << "\nby\n" << efasta(answer) << "\n";
                    PRINT3_TO( rout, TIG, estart, estop );
                    replacements.push( estart, estop, efasta(answer) );    }    }
          return;    }

     // Handle the acyclic case.  First find all paths.

     rout << "acyclic" << endl;
     vec< vec<int> > paths;
     if ( !h.EdgePaths( left_id, right_id, paths, -1, heur.max_paths ) )
          rout << "path count too large\n";
     else
     {    PickPath( low, high, etiglet, tig, estart, estop, h, kbb, basesy, qualsy, 
               Max(readlengths), paths, heur.max_expand_to, heur.q_junk_max, 
               heur.max_non_losers, heur.max_paths, TIG, logc.QLT, tmp_dir, 
               logc.SHOW_VOTES, logc.data_dir, rout, replacements );    }    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_String_OrDefault(SUBDIR, "test");
     CommandArgument_String(SCAFFOLDS_IN);
     CommandArgument_String_OrDefault(SCAFFOLDS_OUT, SCAFFOLDS_IN + ".local");
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
	  "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_Int_OrDefault_Doc(VERBOSITY, 0, "can be 0, 1, or 2");
     CommandArgument_Bool_OrDefault(QLT, False);
     CommandArgument_String_OrDefault_Doc(TIGS, "",
          "if specified, process just these contigs, turns off writing unless "
          "FORCE_WRITE=True");
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(FORCE_WRITE, False);
     CommandArgument_Int_OrDefault(STARTX, -1);
     CommandArgument_Int_OrDefault(STOPX, -1);
     CommandArgument_Bool_OrDefault_Doc(DIRECT, False,
          "log directly to cout, forces VERBOSITY >= 1 and NUM_THREADS = 1");
     CommandArgument_Bool_OrDefault(SHOW_VOTES, False);
     CommandArgument_Bool_OrDefault(DO_CYCLIC, False);
     CommandArgument_Int_OrDefault(WINDOW_START, -1);
     CommandArgument_String_OrDefault_Doc(DUMP_EC, "",
          "file to dump error-corrected reads to, only recommended if used with "
          "TIGS+WINDOW_START so there is only one window analyzed");
     CommandArgument_String_OrDefault_Doc(DOT, "",
          "if specified, name of file to dump DOT graph to; should only be\n"
          "run with a single window");
     CommandArgument_Bool_OrDefault_Doc(RECYCLE, False,
          "reuse reordered read files and don't delete them");
     EndCommandArguments;

     // Define directories, etc.

     double clock = WallClockTime( );
     String data_dir = PRE + "/" + DATA;
     String run_dir = PRE + "/" + DATA + "/" + RUN;
     String sub_dir = run_dir + "/ASSEMBLIES/" + SUBDIR;
     String head = sub_dir + "/" + SCAFFOLDS_IN;
     vec<int> tuse;
     ParseIntSet( TIGS, tuse );
     if ( STARTX >= 0 || STOPX >= 0 ) 
     {    ForceAssertLt( STARTX, STOPX );
          ForceAssert( tuse.solo( ) );    }
     if ( WINDOW_START >= 0 )
     {    ForceAssert( tuse.solo( ) );
          const int default_wing = 6000;
          if ( STARTX < 0 ) STARTX = WINDOW_START - default_wing;
          if ( STOPX < 0 ) STOPX = WINDOW_START + default_wing;    }
     // if ( tuse.solo( ) ) DIRECT = True;

     // Thread control.

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Define heuristics.

     const int K = 20;
     heuristics heur;
     heur.max_reads = 50000;
     heur.min_calls = 10;
     heur.agree_ceil_weak = 0.55;
     heur.agree_floor_strong = 0.9;
     heur.MIN_FRAC = 0.1;
     heur.q_junk_max = 2;
     heur.max_tandem_period = 8;
     heur.min_tandem_copies = 3;
     heur.min_tandem_length = 12;
     heur.max_offby = 10.0;
     heur.max_non_losers = 3;
     heur.max_expand_to = 1000;
     heur.max_paths = 200;
     heur.gapdev_mult = 3.0;
     heur.boundary_push = 20;
     heur.max_cyclic_paths = 5;
     heur.max_cyclic_loops = 100;
     heur.min_devs_off = 10.0;
     heur.do_cyclic = DO_CYCLIC;
     heur.max_cyclic_edges = 50;

     // Define logging control.

     logging_control logc;
     logc.SHOW_VOTES = SHOW_VOTES;
     logc.DUMP_EC = DUMP_EC;
     logc.QLT = QLT;
     logc.data_dir = data_dir;
     logc.DOT = DOT;

     // Load assembly.

     String efasta_file = sub_dir + "/" + SCAFFOLDS_IN + ".contigs.efasta";
     vec<efasta> tigse;
     LoadEfastaIntoStrings( efasta_file, tigse );
     vecbasevector tigs( sub_dir + "/" + SCAFFOLDS_IN + ".contigs.fastb" );
     if ( tuse.empty( ) ) tuse = vec<int>( tigs.size( ), vec<int>::IDENTITY );
     else if ( tuse.back( ) >= (int) tigs.size( ) )
     {    cout << "Illegal value for TIGS, there are only " << tigs.size( )
               << " contigs." << endl;
          exit(1);    }
     String supers_file = sub_dir + "/" + SCAFFOLDS_IN + ".superb";
     vec<superb> scaffolds;
     ReadSuperbs( supers_file, scaffolds );
     const int ploidy = FirstLineOfFile( run_dir + "/ploidy" ).Int( );

     // Build a map of kmers in the contigs.

     cout << Date( ) << ": building kmer maps for contigs" << endl;
     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup( tigs, kmers_plus );
     cout << Date( ) << ": copying" << endl;
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Count number of libraries.

     int nlibs_frag, nlibs_jump, nlibs_long;
     vec<String> lib_types_to_use;
     lib_types_to_use.push_back( "jump" );
     const Bool TIME_STAMPS = True;
     GetLibCounts( run_dir, lib_types_to_use, VERBOSITY, TIME_STAMPS,
		   "","jump_reads_ec","", nlibs_frag, nlibs_jump, nlibs_long );
     int nlibs = nlibs_frag + nlibs_jump + nlibs_long;

     // Align reads.

     cout << Date( ) << ": loading reads" << endl;
     PairsManager pairs( run_dir + "/jump_reads_ec.pairs" );
     vec<unsigned short> read_lengths;
     vec<Bool> placed_fw, placed_rc;
     vec< pair<int,int> > placement;
     vec< vec<longlong> > aligns_index;
     {    vecbasevector jbases( run_dir + "/jump_reads_ec.fastb" );
          read_lengths.resize( jbases.size( ) );
          for ( size_t i = 0; i < jbases.size( ); i++ )
               read_lengths[i] = jbases[i].size( );
          cout << Date( ) << ": aligning reads" << endl;
          const Bool TIME_STAMPS = True;
          AlignReads( 40, 20, True, tigs, jbases, pairs, placed_fw, placed_rc, 
               placement, aligns_index, TIME_STAMPS );    }

     // Define distributions for libraries.

     cout << Date( ) << ": begin defining distributions" << endl;
     vec<Bool> lfail;
     vec<int> max_dist;
     vec< vec<int> > DISTS;
     vec< vec<double> > X;
     const Bool HISTOGRAMS = False;
     DefineDistributions( nlibs_jump, 0, tigs, read_lengths, pairs, placed_fw,
          placed_rc, placement, aligns_index, TIME_STAMPS, HISTOGRAMS, lfail, 
          max_dist, DISTS, X );

     // Set up to access read locations, then build reordered read files.

     read_locs_on_disk locs_file( head, run_dir );
     vec<int64_t> new_frag_id, new_jump_id;
     for ( int pass = 1; pass <= 2; pass++ )
     {    String rhead = ( pass == 1 ? "frag_reads_filt" : "jump_reads_filt" );
          vec<int64_t>& new_id = ( pass == 1 ? new_frag_id : new_jump_id );
          if ( !RECYCLE || !IsRegularFile( head + "." + rhead + ".reorder.ids" ) )
          {    cout << Date( ) << ": loading "
                    << ( pass == 1 ? "frag" : "jump" ) << " read base files" << endl;
               vecbasevector bases( run_dir + "/" + rhead + ".fastb" ), new_bases;
               cout << Date( ) << ": building reordered read base files" << endl;
               new_id.resize( bases.size( ), -1 );
               for ( int tis = 0; tis < tuse.isize( ); tis++ )
               {    vec<read_loc> locs;
                    locs_file.LoadContig( tuse[tis], locs );
                    for ( int j = 0; j < locs.isize( ); j++ )
                    {    const read_loc& rl = locs[j];
                         if ( STARTX >= 0 && rl.Stop( ) <= STARTX ) continue;
                         if ( STOPX >= 0 && rl.Start( ) >= STOPX ) continue;
                         if ( pass == 1 && !rl.Frag( ) ) continue;
                         if ( pass == 2 && !rl.Jump( ) ) continue;
                         new_id[ rl.ReadId( ) ] = new_bases.size( );
                         new_bases.push_back_reserve( 
                              bases[ rl.ReadId( ) ] );    }    }
               new_bases.WriteAll( head + "." + rhead + ".reorder.fastb" );
               BinaryWrite3( head + "." + rhead + ".reorder.ids", new_id );     }
          else BinaryRead3( head + "." + rhead + ".reorder.ids", new_id );
          if ( !RECYCLE || !IsRegularFile( head + "." + rhead + ".reorder.qualb" ) )
          {    cout << Date( ) << ": loading "
                    << ( pass == 1 ? "frag" : "jump" ) << " read qual files" << endl;
               vecqualvector quals( run_dir + "/" + rhead + ".qualb" ), new_quals;
               cout << Date( ) << ": building reordered read qual files" << endl;
               for ( int tis = 0; tis < tuse.isize( ); tis++ )
               {    vec<read_loc> locs;
                    locs_file.LoadContig( tuse[tis], locs );
                    for ( int j = 0; j < locs.isize( ); j++ )
                    {    const read_loc& rl = locs[j];
                         if ( STARTX >= 0 && rl.Stop( ) <= STARTX ) continue;
                         if ( STOPX >= 0 && rl.Start( ) >= STOPX ) continue;
                         if ( pass == 1 && !rl.Frag( ) ) continue;
                         if ( pass == 2 && !rl.Jump( ) ) continue;
                         new_quals.push_back_reserve( 
                              quals[ rl.ReadId( ) ] );    }    }
               new_quals.WriteAll( head + "." + rhead + ".reorder.qualb" );    }    }

     // Go through the contigs.

     vec< vec< pair<int,int> > > markups( tuse.size( ) );
     vec<Bool> DONE( tuse.size( ), False );
     vec<String> REPORT( tuse.size( ) );
     int REPORT_INDEX = 0;
     unsigned tigs_processed = 0;
     cout << Date( ) << ": processing " << tuse.size( ) << " contigs" << endl;
     if (DIRECT) omp_set_num_threads(1);
     #pragma omp parallel for schedule(dynamic, 1)
     for ( int tis = 0; tis < tuse.isize( ); tis++ )
     {    int TIG = tuse[tis];
          if ( VERBOSITY == 0 && tuse.size( ) > 1 ) {
               #pragma omp critical
               dots_pct(tigs_processed++, tuse.size());
          }

          basevector tig = tigs[TIG];
          ostringstream toutx;
          ostream& tout = ( DIRECT ? cout : toutx );
          tout << "\n" << Date( ) << ": analyzing contig " << TIG << " of " 
               << tuse.size( ) << ", l = " << tig.size( ) << " bases" << endl;
          int START = 0, STOP = tig.size( );
          if ( STARTX >= 0 ) START = STARTX;
          if ( STOPX >= 0 ) STOP = STOPX;

          // Load read locations and remove those that are not within bounds.

          vec<read_loc> locs;
          tout << Date( ) << ": loading read locations" << endl;
          #pragma omp critical
          {    locs_file.LoadContig( TIG, locs );    }
          vec<Bool> locs_delete( locs.size( ), False );
          for ( int j = 0; j < locs.isize( ); j++ )
          {    if ( START >= 0 && locs[j].Stop( ) <= START ) locs_delete[j] = True;
               if ( STOP >= 0 && locs[j].Start( ) >= STOP ) 
                    locs_delete[j] = True;    }
          EraseIf( locs, locs_delete );
     
          // Sort out read locations.

          vec<int> ids_frag, ids_jump;
          vec<Bool> or_frag, or_jump;
          vec<int> loc_id_frag, loc_id_jump;
          for ( int i = 0; i < locs.isize( ); i++ )
          {    const read_loc& rl = locs[i];
               if ( rl.Frag( ) ) 
               {    ids_frag.push_back( new_frag_id[ rl.ReadId( ) ] );
                    or_frag.push_back( rl.Fw( ) );
                    loc_id_frag.push_back(i);    }
               if ( rl.Jump( ) ) 
               {    ids_jump.push_back( new_jump_id[ rl.ReadId( ) ] );
                    or_jump.push_back( rl.Fw( ) );
                    loc_id_jump.push_back(i);    }    }
          SortSync( ids_frag, loc_id_frag );
          SortSync( ids_jump, loc_id_jump );

          // Load reads.

          tout << Date( ) << ": loading reads" << endl;
          vecbasevector bases;
          bases.Read( head + ".frag_reads_filt.reorder.fastb", ids_frag );
          bases.Read( head + ".jump_reads_filt.reorder.fastb", ids_jump );
          vecqualvector quals;
          quals.Read( head + ".frag_reads_filt.reorder.qualb", ids_frag );
          quals.Read( head + ".jump_reads_filt.reorder.qualb", ids_jump );
     
          // Put reads in contig orientation.
     
          tout << Date( ) << ": orienting pairs" << endl;
          for ( int i = 0; i < ids_frag.isize( ); i++ )
          {    if ( !locs[ loc_id_frag[i] ].Fw( ) )
               {    bases[i].ReverseComplement( );
                    quals[i].ReverseMe( );    }    }
          for ( int i = 0; i < ids_jump.isize( ); i++ )
          {    if ( !locs[ loc_id_jump[i] ].Fw( ) )
               {    bases[ i + ids_frag.isize( ) ].ReverseComplement( );
                    quals[ i + ids_frag.isize( ) ].ReverseMe( );    }    }
     
          // Create data structure to store signal.
          // Current criterion to flag base as bad, one or more of the following:
          // (a) pileup of >= 10 bases at < 50% concordance with reference
          // (b) ambiguity
          // (c) repeat (not implemented).

          tout << Date( ) << ": looking for signal" << endl;
          vec< triple<int,int,String> > signal;
          vec<Bool> weak;

          // Scan to locate ambiguities.

          int pos = 0;
          efasta& T = tigse[TIG];
          for ( int i = 0; i < T.isize( ); i++ )
          {    if ( T[i] == '{' )
               {    int start = pos, stop = pos, istart = i;
                    for ( i++; i < T.isize( ); i++ )
                    {    if ( T[i] == ',' ) break;
                         stop++; pos++;    }
                    for ( i++; i < T.isize( ); i++ )
                         if ( T[i] == '}' ) break;
                    if ( IntervalOverlap( start, stop, START, STOP ) > 0 )
                    {    signal.push( start, stop, "ambiguity "
                              + T.substr( istart, i + 1 - istart ) );    
                         weak.push_back(True);    }    }
               else pos++;    }

          // Scan to locate repeats.
     
          for ( int cpos = 0; cpos < tig.isize( ); cpos++ )
          {    for ( int period = 1; period <= heur.max_tandem_period; period++ )
               {    if ( cpos + heur.min_tandem_copies * period > tig.isize( ) ) 
                         break;
                    int len;
                    for ( len = 1; cpos + len < tig.isize( ); len++ )
                         if ( tig[cpos + len] != tig[cpos + len % period] ) break;
                    if ( len < heur.min_tandem_length ) continue;
                    if ( len < heur.min_tandem_copies * period ) continue;
                    signal.push( cpos, cpos + len, 
                         "tandem repeat, period " + ToString(period) );
                    weak.push_back(True); 
                    cpos += len - 1;    }    }
      
          // Form pileup, identify weak positions and strong windows.

          vec<dumbcall> calls, qcalls;
          calls.resize_and_set( tig.size( ), dumbcall( ) );
          qcalls.resize_and_set( tig.size( ), dumbcall( ) );
          for ( int i = 0; i < ids_frag.isize( ); i++ )
          {    const basevector& b = bases[i];
               const qualvector& q = quals[i];
               AddToPileup( locs[ loc_id_frag[i] ], b, tig, calls );
               AddToPileup( locs[ loc_id_frag[i] ], b, q, tig, qcalls );    }
          for ( int i = 0; i < ids_jump.isize( ); i++ )
          {    const basevector& b = bases[ i + ids_frag.isize( ) ];
               const qualvector& q = quals[ i + ids_frag.isize( ) ];
               AddToPileup( locs[ loc_id_jump[i] ], b, tig, calls );
               AddToPileup( locs[ loc_id_jump[i] ], b, q, tig, qcalls );    }
          for ( int j = START; j < STOP; j++ )
          {    int ncalls = 0, nqcalls = 0, nqcalls_agree = 0;
               for ( int k = 0; k < 6; k++ )
               {    ncalls += calls[j].base[k];
                    nqcalls += qcalls[j].base[k];
                    if ( k == tig[j] ) nqcalls_agree += qcalls[j].base[k];    }
               double qagree = double(nqcalls_agree)/double(nqcalls);
               if ( ncalls >= heur.min_calls && qagree < heur.agree_ceil_weak ) 
               {    signal.push( 
                         j, j+1, "agree = " + ToString( 100.0 * qagree ) + "%" );
                    weak.push_back(True);    }    }
          vec<Bool> marked( STOP - START, False );
          for ( int i = 0; i < signal.isize( ); i++ )
          {    for ( int j = signal[i].first; j < signal[i].second; j++ )
                    if ( j >= START && j < STOP ) marked[ j - START ] = True;    }
          vec<Bool> strong( STOP - START, False );
          for ( int j = START; j < STOP; j++ )
          {    int ncalls = 0, nqcalls = 0, nqcalls_agree = 0;
               for ( int k = 0; k < 6; k++ )
               {    ncalls += calls[j].base[k];
                    nqcalls += qcalls[j].base[k];
                    if ( k == tig[j] ) nqcalls_agree += qcalls[j].base[k];    }
               double qagree = double(nqcalls_agree)/double(nqcalls);
               if ( ncalls >= heur.min_calls && qagree >= heur.agree_floor_strong
                    && !marked[ j - START ] )
               {    strong[ j - START ] = True;    }    }
          for ( int i = 0; i < strong.isize( ); i++ )
          {    if ( !strong[i] ) continue;
               int j;
               for ( j = i + 1; j < strong.isize( ); j++ )
                    if ( !strong[j] ) break;
               if ( j - i >= K ) 
               {    signal.push( START + i, START + j, "" );
                    weak.push_back(False);    }
               i = j - 1;    }

          // Define windows.  Note that windows are in fasta/fastb coordinates, not 
          // in efasta coordinates.

          SortSync( signal, weak );
          vec< pair<int,int> > windows;
          vec<String> reports;
          for ( int i = 0; i < signal.isize( ); i++ )
          {    if ( weak[i] )
               {    int j;
                    for ( j = i + 1; j < signal.isize( ); j++ )
                         if ( !weak[j] ) break;
                    if ( i > 0 && j < signal.isize( ) )
                    {    int start = signal[i-1].second - K;
                         int stop = signal[j].first + K;
                         windows.push( start, stop );
                         ostringstream out;
                         for ( int l = i; l < j; l++ )
                         {    out << TIG << "." << signal[l].first << "-" 
                                   << signal[l].second << "   " << signal[l].third 
                                   << "\n";    }
                         reports.push_back( out.str( ) );    }
                    i = j - 1;    }    }

          // Get maximum read length.

          vec<int> readlengths;
          for ( size_t i = 0; i < bases.size( ); i++ )
               readlengths.push_back( bases[i].size( ) );
          UniqueSort(readlengths);
          tout << Date( ) << ": nreadlengths = " << readlengths.size( ) << endl;

          // Go through the windows.

          vec<String> report( windows.size( ) );
          vec< triple<int,int,int> > markupsx( windows.size( ) );
          vec<Bool> markupsx_used( windows.size( ), False );
          vec< triple<int,int,String> > replacements;
          vec< triple<int,int,String> > replacementsx( windows.size( ) );
          vec<Bool> replacementsx_used( windows.size( ), False );
          vec<Bool> done( windows.size( ), False );
          tout << Date( ) << ": traversing " << windows.size( ) << " windows" 
               << endl;
          double clock = WallClockTime( );
          redo:
          for ( int i = 0; i < windows.isize( ); i++ )
          {    if ( done[i] ) continue;
               if ( WINDOW_START >= 0 && windows[i].first != WINDOW_START ) continue;
               if (DIRECT)
               {    cout << "\n=============================================="
                         << "======================================\n";
                    cout << "\n" << TIG << "." << windows[i].first << "-" 
                         << windows[i].second << "   WINDOW " << i << "\n";    }
               ostringstream routx;
               ostream& rout = ( DIRECT ? cout : routx );
               vec< triple<int,int,int> > markups_this;
               vec< triple<int,int,String> > replacements_this;
               double wclock = WallClockTime( );
               ProcessWindow( kmers_plus, kmers, heur, logc, 
                    head + ".FixLocal." + ToString(i) + "x",
                    windows[i], reports[i], tig, T, TIG, K, ids_frag, 
                    ids_jump, locs, loc_id_frag, loc_id_jump, bases, quals, 
                    readlengths, tigs, aligns_index, placed_fw, placed_rc, 
                    placement, nlibs_jump, 0, DISTS, X, lfail, read_lengths, 
                    pairs, rout, markups_this, replacements_this );
               rout << "\ntime used for window = " << TimeSince(wclock) << endl;
               if ( !DIRECT ) report[i] = routx.str( );
               if ( markups_this.nonempty( ) ) 
               {    markupsx[i] = markups_this[0];
                    markupsx_used[i] = True;    }
               if ( replacements_this.nonempty( ) )
               {    replacementsx[i] = replacements_this[0];    
                    replacementsx_used[i] = True;    }    }

          // Print reports.

          if ( VERBOSITY >= 2 && !DIRECT )
          {    for ( int i = 0; i < report.isize( ); i++ )
               {    if ( !done[i] ) 
                    {    tout << "\n=============================================="
                              << "======================================\n";
                         tout << "\n" << TIG << "." << windows[i].first << "-" 
                              << windows[i].second << "   WINDOW " << i << "\n";
                         tout << report[i];    }    }    }
          for ( int i = 0; i < report.isize( ); i++ )
               done[i] = True;

          // Check for overlapping windows.  Note that some of this is unnecessary:
          // the original windows overlap in some cases, which didn't have to be
          // the case.  On the other hand, by pushing the boundaries of some 
          // windows, we will occasionally get overlaps that could not have been
          // prevented.

          vec<Bool> to_remove( windows.size( ), False );
          for ( int i = 0; i < windows.isize( ) - 1; i++ )
          {    if ( to_remove[i] ) continue;
               for ( int j = i+1; j < windows.isize( ); j++ )
               {    if ( to_remove[j] ) continue;
                    if ( ( markupsx_used[i] || replacementsx_used[i] )
                         && ( markupsx_used[j] || replacementsx_used[j] )
                         && windows[i].second >= windows[j].first )
                    {    int low = Min( windows[i].first, windows[j].first );
                         int high = Max( windows[i].second, windows[j].second );
                         windows[i].first = low, windows[i].second = high;
                         reports[i] = reports[i] + reports[j];
                         done[i] = False;
                         to_remove[j] = True;    }    }    }
          int redos = done.isize( ) - Sum(done);
          if ( redos > 0 )
          {    if ( VERBOSITY >= 2 )
               {    tout << "\n" << Date( ) << ": detected " << redos << " overlaps"
                         << " between windows, redoing some" << endl;    }
               EraseIf( windows, to_remove );
               EraseIf( done, to_remove );
               EraseIf( report, to_remove );
               EraseIf( reports, to_remove );
               EraseIf( markupsx, to_remove );
               EraseIf( markupsx_used, to_remove );
               EraseIf( replacementsx, to_remove );
               EraseIf( replacementsx_used, to_remove );
               goto redo;    }
          tout << Date( ) << ": done traversing windows, time used = " 
               << TimeSince(clock) << endl;

          // Push back results.

          for ( int j = 0; j < windows.isize( ); j++ )
          {    if ( replacementsx_used[j] ) 
                    replacements.push_back( replacementsx[j] );
               if ( markupsx_used[j] )
                    markups[tis].push( markupsx[j].second, markupsx[j].third );    }
          tout << Date( ) << ": made " << replacements.size( ) << " replacements "
               << "and " << markups[tis].size( ) << " markups" << endl;
	  
          // Insert replacements.
	  
          ReverseSort(replacements);
          for ( int i = 0; i < replacements.isize( ); i++ ){    
	    int estart = replacements[i].first, estop = replacements[i].second;
	    String left = T.substr( 0, estart );
	    String middle = replacements[i].third;
	    String right = T.substr( estop, T.isize( ) - estop );

	    T = left + middle + right;
            
	    ValidateEfastaRecord( T, 
				  "in processing of contig " + ToString(TIG) );    
	    
	    int deltaoff = middle.size() - (estop - estart);
	    for ( int j = 0; j < markups[tis].isize(); j++ ){
	      int& m0 = markups[tis][j].first;
	      int& m1 = markups[tis][j].second;
              #pragma omp critical
	      { 
		if ( m0 >= estop && m1 >= estop){
		  m0 += deltaoff;
		  m1 += deltaoff;
		}else if ( m0 >= estart && m1 >= estop ){
		  m0  = estop + deltaoff; // shorten
		  m1 += deltaoff;
		}
		else if ( m0 >= estart && m1 < estop ){
		  m1 = m0; // remove
		}
		else if ( m0 < estart && m1 >= estart && m1 < estop){
		  m1 = estart -1; // shorten
		}
		else if ( m0 < estart && m1 >= estop){
		  m1 += deltaoff;
		}
	      }
	    }
	  }
	  
	  // Dump report.
	  
	  DONE[tis] = True, REPORT[tis] = toutx.str( );
	  if ( VERBOSITY >= 1 )
	    {
            #pragma omp critical
	      {    for ( ; REPORT_INDEX < (int) tuse.size( ); REPORT_INDEX++ )
		  {    if ( !DONE[REPORT_INDEX] ) break;
		    cout << REPORT[REPORT_INDEX];    
		    flush(cout);    }    }    }    
     }
     
     // Clean up and write new assembly.
     cout << "Done with fixing, preparing output" << endl;

     if ( !RECYCLE )
     {    for ( int pass = 1; pass <= 2; pass++ )
          {    String rhead = ( pass == 1 ? "frag_reads_filt" : "jump_reads_filt" );
               Remove( head + "." + rhead + ".reorder.fastb" );
               Remove( head + "." + rhead + ".reorder.ids" );
               Remove( head + "." + rhead + ".reorder.qualb" );    }    }
     if ( WRITE && ( TIGS == "" || FORCE_WRITE ) )
     {    if ( VERBOSITY > 0 ) cout << "\n";
          cout << Date( ) << ": writing assembly" << endl;
          for ( int s = 0; s < scaffolds.isize( ); s++ )
          {    superb& S = scaffolds[s];
               for ( int j = 0; j < S.Ntigs( ); j++ )
                    S.SetLen( j, tigse[ S.Tig(j) ].Length1( ) );    }
          Assembly A( scaffolds, tigse );
          A.WriteAll( sub_dir + "/" + SCAFFOLDS_OUT );
          Ofstream( out, sub_dir + "/" + SCAFFOLDS_OUT + ".markup" );
          for ( int tis = 0; tis < tuse.isize( ); tis++ ){    
	    for ( int j = 0; j < markups[tis].isize( ); j++ ){    
	      if ( markups[tis][j].second > markups[tis][j].first && 
		   markups[tis][j].first >= 0 &&
		   markups[tis][j].second < tigse[ tuse[tis] ].isize() )
		out << tuse[tis] << " " << markups[tis][j].first
		    << " " << markups[tis][j].second << "\n";    
	    }    
	  }    
     }
     cout << "\n" << Date( ) << ": done, total time used = " << TimeSince(clock) 
          << endl;    }
