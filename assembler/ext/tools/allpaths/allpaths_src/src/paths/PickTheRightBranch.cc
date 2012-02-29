///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// PickTheRightBranch.  At junctions in the unipath graph, pick the right branch, 
// then murder the other one.

// MakeDepend: dependency QueryLookupTable

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "PairsManager.h"
#include "Qualvector.h"
#include "Set.h"
#include "VecUtilities.h"
#include "graph/Digraph.h"
#include "kmers/KmerParcelsTools.h"
#include "math/Functions.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/GetNexts.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/PdfEntry.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/UnibaseCopyNumberCommon.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"
#include "paths/UnipathFixerTools.h"
#include "random/Bernoulli.h"
#include "util/SearchFastb2Core.h"

int main( int argc, char** argv ) 
{
  
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String( PRE );
     CommandArgument_String( DATA );
     CommandArgument_String( RUN );
     CommandArgument_Int( K );
     CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
       "Number of threads to use (use all available processors if set to 0)");
     CommandArgument_String_OrDefault( UNIBASES_IN, "all_reads" );
     CommandArgument_String_OrDefault( READS_EC_IN, "frag_reads_edit" );
     CommandArgument_String_OrDefault( UNIBASES_OUT, "all_reads.picked" );
     CommandArgument_Int_OrDefault( MIN_KMERS, 0 );
     CommandArgument_Bool_OrDefault( VALIDATE, False );
     CommandArgument_Bool_OrDefault( STATS, False );
     CommandArgument_Bool_OrDefault( PRINT_READ_IDS, False );
     CommandArgument_Bool_OrDefault( CN_REACH, True );
     CommandArgument_Bool_OrDefault( CN_ALIGN, False );
     CommandArgument_Bool_OrDefault( WRITE, True );
     CommandArgument_Int_OrDefault( VERBOSITY, 0 );
     EndCommandArguments;

     // Thread control
     
     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );
  
     // Define heuristics.

     const double max_p_opp = 0.001;

     // Define directories.
  
     String data_dir = PRE + "/" + DATA;
     String run_dir = data_dir + "/" + RUN;

     // Load unibases and ancillary data.

     String unibases_head = run_dir + "/" + UNIBASES_IN + ".unibases";
     String kK = ".k" + ToString(K);
     vecbasevector unibases( unibases_head + kK );
     int nuni = unibases.size( );
     // VecPdfEntryVec CNs(run_dir + "/" + READS + ".unipaths.predicted_count" + kK);
     vec<int> predicted_CNs( nuni, -1 ), to_rc;
     // for ( size_t i = 0; i < nuni; i++ )
     //      GetMostLikelyValue( predicted_CNs[i], CNs[i] );
     UnibaseInvolution( unibases, to_rc, K );

     // Determine which unibases are genomic.

     vec<int> true_CN(nuni, 0);
     if (VALIDATE)
     {    vec< triple<int64_t,int64_t,int> > ALIGNS;
          SearchFastb2( unibases_head + kK, data_dir + "/genome.fastb", K,
                          &ALIGNS, 0, -1, 0.9, False );
          for ( size_t i = 0; i < ALIGNS.size( ); i++ )
               true_CN[ ALIGNS[i].first ]++;    }

     // Define the graph U whose vertices are unibases.

     digraph U;
     {    vec< vec<int> > nexts, backs(nuni);
          GetNexts( K, unibases, nexts );
          for ( int u = 0; u < nuni; u++ )
          {    Sort( nexts[u] );
               for ( int j = 0; j < nexts[ to_rc[u] ].isize( ); j++ )
                    backs[u].push_back( to_rc[ nexts[ to_rc[u] ][j] ] );
               Sort( backs[u] );    }
          U.Initialize( nexts, backs );    }
     HyperBasevector h;
     BuildUnibaseAdjacencyHyperBasevector( K, U, unibases, h );
     vec<int> L;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          L.push_back( h.EdgeObject(i).size( ) - K + 1 );
     digraphE<int> G( h, L );

     // Load read pairs and quality scores.

     String reads_head = run_dir + "/" + READS_EC_IN;
     PairsManager pairs( reads_head + ".pairs" );
     vecqualvector quals( reads_head + ".qualb" );
     int64_t nfrags = quals.size( );

     vecbasevector bases0( run_dir + "/frag_reads_filt.fastb" );
     vecqualvector quals0( run_dir + "/frag_reads_filt.qualb" );
     vecbasevector jbases0( run_dir + "/jump_reads_filt.fastb" );
     vecqualvector jquals0( run_dir + "/jump_reads_filt.qualb" );

     // Build paths database for the reads and unibases.  Note:
     // (1) Locating the database creation here is temporary.
     // (2) It may already exist.
     // (3) It is not completely clear which reads should be used.  For now we use
     //     frag_reads_edit.
     // (4) We could use jump reads too but if they have been extended they might
     //     pollute the counts.

     vecbasevector all(unibases);
     all.ReadAll( reads_head + ".fastb", True );
     all.ReadAll( run_dir + "/jump_reads_filt.fastb", True );
     vecKmerPath paths, pathsrc;
     vec<tagged_rpint> pathsdb;
     ReadsToPathsCoreY( all, K, paths, pathsrc, pathsdb,
          run_dir + "/PickTheRightBranch", NUM_THREADS );
     vecKmerPath unipaths(paths);
     unipaths.resize(nuni);
     vec<tagged_rpint> unipathsdb;
     CreateDatabase( unipaths, unipathsdb );

     // Find the unipaths that live between the ends of fragment pairs, and create
     // "extra_hits" corresponding to them.  Note that this loop could be 
     // parallelized, although it appears to be very fast.

     cout << Date( ) << ": looking between fragment pairs" << endl;
     vec<int> extra_hits( nuni, 0 );
     {    vec<int> to_left, to_right;
          G.ToLeft(to_left), G.ToRight(to_right);
          for ( size_t pid = 0; pid < pairs.nPairs( ); pid++ )
          {    int id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
               if ( id1 < 0 || id2 < 0 ) continue;
               longlong id1_end = paths[nuni+id1].LastSegment( ).Stop( );
               int u1 = -1, u2 = -1;
               vec<longlong> con;
               Contains( unipathsdb, id1_end, con );
               for ( int l = 0; l < con.isize( ); l++ )
               {    const tagged_rpint& t = unipathsdb[ con[l] ];
                    u1 = t.ReadId( );
                    break;    }
               longlong id2_begin = paths[nuni+id2].FirstSegment( ).Start( );
               Contains( unipathsdb, id2_begin, con );
               for ( int l = 0; l < con.isize( ); l++ )
               {    const tagged_rpint& t = unipathsdb[ con[l] ];
                    u2 = to_rc[ t.ReadId( ) ];
                    break;    }
               if ( u1 < 0 || u2 < 0 || u1 == u2 ) continue;
               if ( Member( U.From(u1), u2 ) ) continue;
               if ( VERBOSITY >= 1 )
               {    cout << "\n";
                    PRINT4( u1, u2, id1, id2 );    }
               const int max_paths = 1000;
               int max_path_length = pairs.sep(pid) + 3 * pairs.sd(pid) + K - 1;
               vec< vec<int> > upaths;
               if ( G.AllPathsLengthRange( to_right[u1], to_left[u2], 
                    0, max_path_length, to_right, upaths, max_paths ) )
               {    vec<int> common;
                    if ( upaths.nonempty( ) ) Intersection( upaths, common );
                    for ( int r = 0; r < common.isize( ); r++ )
                    {    extra_hits[ common[r] ]++;
                         if ( to_rc[ common[r] ] != common[r] )
                              extra_hits[ to_rc[ common[r] ] ]++;    }
                    if ( VERBOSITY >= 1 )
                    {    cout << "common:";
                         for ( int r = 0; r < common.isize( ); r++ )
                              cout << " " << common[r];
                         cout << "\n";    }    }    }    }
     cout << "\n" << Date( ) << ": done" << endl;

     // Find the reads that do not share a K-mer with a unibase.  Find the longest
     // perfect alignment to a unibase, that is at least 40 bases long.  If there
     // is exactly one such alignment, record the read as providing coverage for
     // that unibase.  Note that the case where the read lands on the overlap
     // between two unibases is excluded, and perhaps this isn't optimal.  Note
     // also that these extra reads are currently used only to improve copy number
     // computations, and not to supplement counts for picking branches.
     //
     // This code doesn't seem to help much but it's an interesting concept, maybe
     // useful later.

     vec<int> extra( nuni, 0 );
     if (CN_ALIGN)
     {    const int seed = 40;
          String xreads_file = reads_head + ".fastb.uncorrected";
          {    vecbasevector xreads( reads_head + ".fastb" );
               for ( size_t id = 0; id < xreads.size( ); id++ )
               {    xreads[id].resize(seed);
                    for ( int j = 0; j < paths[id].NSegments( ); j++ )
                    {    vec<longlong> con;
                         Contains( pathsdb, paths[id].Segment(j), con );
                         for ( int l = 0; l < con.isize( ); l++ )
                         {    const tagged_rpint& t = pathsdb[ con[l] ];
                              if ( t.ReadId( ) < nuni )
                              {    xreads[id].resize(0);
                                   goto tail;    }    }    }
                    tail: continue;    }
               xreads.WriteAll(xreads_file);    }
          vec< triple<int64_t,int64_t,int> > ALIGNS;
          SearchFastb2( xreads_file, unibases_head + kK, seed, &ALIGNS, 0,
               -1, 0.9, True );
          Remove(xreads_file);
          {    vec<Bool> aligns_to_remove( ALIGNS.size( ), False );
               for ( size_t i = 0; i < ALIGNS.size( ); i++ )
                    if ( ALIGNS[i].third < 0 ) aligns_to_remove[i] = True;
               EraseIf( ALIGNS, aligns_to_remove );    }
          for ( size_t i = 0; i < ALIGNS.size( ); i++ )
          {    size_t j; 
               int64_t id = ALIGNS[i].first;
               for ( j = i + 1; j < ALIGNS.size( ); j++ )
                    if ( ALIGNS[j].first != id ) break;
               vec<int> align_length( j - i, seed );
               vec<size_t> align_ids;
               for ( size_t l = i; l < j; l++ )
                    align_ids.push_back(l);
               for ( size_t l = i; l < j; l++ )
               {    size_t u = ALIGNS[i].second;
                    int pos = ALIGNS[i].third;
                    const basevector &rd1 = all[ nuni + id ], &rd2 = unibases[u];
                    for ( int r = seed; r < K; r++ )
                    {    if ( r >= rd1.isize( ) || pos+r >= rd2.isize( ) ) break;
                         if ( rd1[r] == rd2[pos+r] ) ++align_length[l-i];
                         else break;    }    }
               ReverseSortSync( align_length, align_ids );
               int max_align_length = Max(align_length);
               size_t nkeep;
               for ( nkeep = 0; nkeep < j-i; nkeep++ )
                    if ( align_length[nkeep] < max_align_length ) break;
               align_ids.resize(nkeep);
               if ( nkeep == 1 )
               {    int u = ALIGNS[ align_ids[0] ].second;
                    extra[u]++;
                    extra[ to_rc[u] ]++;    }
               i = j - 1;    }    }

     // Directly compute copy number.  Note that this computation is NOT right.
     // Since we count the number of reads incident upon a given unipath, the 
     // counts are exaggerated for very small unipaths.  This could be fixed with
     // a correction factor.

     int PLOIDY = StringOfFile( data_dir + "/ploidy", 1 ).Int( );
     int64_t total_ids = 0, total_kmers = 0;
     double CN_error = 0.0;
     int CN_error_count = 0;
     vec<double> CN_float(nuni);
     vec<int> count(nuni);
     for ( int pass = 1; pass <= 2; pass++ )
     {    int min_length, nlongest;
          if ( pass == 1 ) GetMinLength( unibases, min_length, nlongest );
          for ( int u = 0; u < nuni; u++ )
          {    if ( pass == 1 && unibases[u].isize( ) < min_length ) continue;

               // Define an extended version of u.  Suppose that the unipaths 
               // x1...xr that go into u do not go to any other unipaths.  Let m_x
               // be the minimum number of kmers in x1...xr.  Then the extended
               // version of u includes the last m_x kmers in each of x1...xr.
               // Similarly, extend on the other side if possible.

               vec< triple<int,int,int> > uext;
               uext.push( u, 0, paths[u].KmerCount( ) );
               int uext_kmers = paths[u].KmerCount( );
               if (CN_REACH)
               {    Bool alt_out = False, alt_in = False;
                    int min_in = 1000000000, min_out = 1000000000;
                    for ( int j = 0; j < U.To(u).isize( ); j++ )
                    {    int x = U.To(u)[j];
                         if ( x == u || !U.From(x).solo( ) ) alt_out = True;
                         min_in = Min( min_in, paths[x].KmerCount( ) );    }
                    if ( !alt_out && min_in < 1000000000 )
                    {    for ( int j = 0; j < U.To(u).isize( ); j++ )
                         {    int x = U.To(u)[j];
                              int pk = paths[x].KmerCount( );
                              uext.push( x, pk - min_in, pk );    }
                         uext_kmers += min_in;    }
                    for ( int j = 0; j < U.From(u).isize( ); j++ )
                    {    int y = U.From(u)[j];
                         if ( y == u || !U.To(y).solo( ) ) alt_in = True;
                         min_out = Min( min_out, paths[y].KmerCount( ) );    }
                    if ( !alt_in && min_out < 1000000000 )
                    {    for ( int j = 0; j < U.From(u).isize( ); j++ )
                         {    int y = U.From(u)[j];
                              int pk = paths[y].KmerCount( );
                              uext.push( y, 0, min_out );    }
                         uext_kmers += min_out;    }    }
               cout << "u copy number, using:\n"; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
               for ( int z = 0; z < uext.isize( ); z++ ) // XXXXXXXXXXXXXXXXXXXXXXXX
               {    int v = uext[z].first; // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
                    int vstart = uext[z].second, vstop = uext[z].third; // XXXXXXXXX
                    PRINT3( v, vstart, vstop );    } // XXXXXXXXXXXXXXXXXXXXXXXXXXXX

               // Find reads incident upon the extended version of u.

               vec<int> ids;
               for ( int z = 0; z < uext.isize( ); z++ )
               {    int v = uext[z].first;
                    int vstart = uext[z].second, vstop = uext[z].third;
                    int pos = 0;
                    for ( int s = 0; s < paths[v].NSegments( ); s++ )
                    {    longlong start = paths[v].Start(s), stop = paths[v].Stop(s);
                         int pos0 = pos;
                         pos += stop - start + 1;
                         if ( pos0 < vstart ) start += vstart - pos0;
                         if ( pos0 + (stop-start+1) > vstop )
                              stop -= ( pos0 + (stop-start+1) - vstop );
                         if ( start <= stop )
                         {    vec<longlong> con;
                              Contains( pathsdb, KmerPathInterval(start,stop), con );
                              for ( int l = 0; l < con.isize( ); l++ )
                              {    const tagged_rpint& t = pathsdb[ con[l] ];
                                   int id = t.ReadId( );
                                   if ( id >= (int) nuni ) 
                                        ids.push_back( id-nuni );    }    }    }    }

               // Compute copy number.

               UniqueSort(ids);
               int nids = ids.size( );
               if (CN_ALIGN) nids += extra[u];
               count[u] = nids;
               if ( pass == 1 )
               {    total_ids += nids, total_kmers += uext_kmers;    }
               else
               {    double CN = double(PLOIDY) 
                         * ( double(nids) / double(uext_kmers) )
                         / ( double(total_ids) / double(total_kmers) );
                    CN_float[u] = CN;
                    int kmers = paths[u].KmerCount( );
                    PRINT5( u, kmers, CN, true_CN[u], ids.size( ) );
                    if ( VALIDATE && true_CN[u] > 0 ) 
                    {    CN_error += Abs( CN - true_CN[u] );
                         CN_error_count++;    }
                    predicted_CNs[u] = int(round(CN));    }    }    }

     // Raise copy numbers.  Suppose that all edges going into u go nowhere else.
     // Compute the number of reads going into u that are supported by at least 20
     // reads.  We raise the copy number of u to at least this number.  Ditto for
     // the other direction.  Note that this is probably too weak in the ploidy 2
     // case.
     //
     // On a 5M test case, this helped slightly, but the effect was small, so I've
     // commented it out.

     /*
     const int min_count_to_qualify = 20;
     {    String unipaths_head = run_dir + "/" + READS + "." + UNIPATHS;
          vecKmerPath unipaths( unipaths_head + kK );
          HyperKmerPath h;
          BuildUnipathAdjacencyHyperKmerPath( K, U, unipaths, h );
          for ( int x = 0; x < h.N( ); x++ )
          {    if ( h.From(x).solo( ) )
               {    int u = h.EdgeObjectIndexByIndexFrom(x, 0);
                    int cn = 0;
                    for ( int j = 0; j < h.To(x).isize( ); j++ )
                    {    int v = h.EdgeObjectIndexByIndexTo(x, j);
                         if ( count[v] >= min_count_to_qualify ) cn++;    }
                    if ( predicted_CNs[u] <= true_CN[u] && cn > true_CN[u] )
                    {    cout << "raising copy number of " << u << " to "
                              << cn << ", true value = " << true_CN[u] << endl;    }
                    predicted_CNs[u] = Max( predicted_CNs[u], cn );    }
               if ( h.To(x).solo( ) )
               {    int u = h.EdgeObjectIndexByIndexTo(x, 0);
                    int cn = 0;
                    for ( int j = 0; j < h.From(x).isize( ); j++ )
                    {    int v = h.EdgeObjectIndexByIndexFrom(x, j);
                         if ( count[v] >= min_count_to_qualify ) cn++;    }
                    if ( predicted_CNs[u] <= true_CN[u] && cn > true_CN[u] )
                    {    cout << "raising copy number of " << u << " to "
                              << cn << ", true value = " << true_CN[u] << endl;    }
                    predicted_CNs[u] = Max( predicted_CNs[u], cn );    }    }    }
     */
      
     // Assess copy numbers.

     if (VALIDATE)
     {    cout << "mean error for floating point copy number: "
               << CN_error / double(CN_error_count) << endl;
          int CN1_wrong = 0, CN1_total = 0, CN_high_wrong = 0, CN_high_total = 0;
          for ( int u = 0; u < nuni; u++ )
          {    if ( true_CN[u] == 1 )
               {    CN1_total++;
                    if ( predicted_CNs[u] > 1 ) CN1_wrong++;    }
               if ( true_CN[u] > 1 )
               {    CN_high_total++;
                    if ( predicted_CNs[u] <= 1 ) CN_high_wrong++;    }    }
          cout << "fraction of true CN 1 reported as CN > 1: "
               << PERCENT_RATIO( 3, CN1_wrong, CN1_total ) << endl;
          cout << "fraction of true CN > 1 reported as CN <= 1: "
               << PERCENT_RATIO( 3, CN_high_wrong, CN_high_total ) << endl;
          for ( int pass = 1; pass <= 2; pass++ )
          {    cout << "\nmedian reported CN ";
               if ( pass == 1 ) cout << "for true CN=1, ";
               cout << "as function of nkmers:\n";
               vec< pair<int,double> > kmers_cn;
               for ( int u = 0; u < nuni; u++ )
               {    if ( true_CN[u] == 1 || pass == 2 )
                         kmers_cn.push( paths[u].KmerCount( ), CN_float[u] );    }
               Sort(kmers_cn);
               const int min_bin = 100;
               for ( int i = 0; i < kmers_cn.isize( ); i++ )
               {    int j, ilast = i;
                    need_more:
                    for ( j = ilast; j < kmers_cn.isize( ); j++ )
                         if ( kmers_cn[j].first != kmers_cn[ilast].first ) break;
                    ilast = j;
                    if ( j - i < min_bin && j < kmers_cn.isize( ) ) goto need_more;
                    if ( j - i < min_bin ) break;
                    vec<double> cn;
                    for ( int l = i; l < j; l++ )
                         cn.push_back( kmers_cn[l].second );
                    Sort(cn);
                    cout << ( kmers_cn[j-1].first + kmers_cn[i].first ) / 2 << " "
                         << cn[ cn.size( ) / 2 ] << "\n";
                    i = j - 1;    }    }    }

     // Compute copy number conservation error at junctions.

     /*
     String unipaths_head = run_dir + "/" + READS + "." + UNIPATHS;
     vecKmerPath unipaths( unipaths_head + kK );
     HyperKmerPath h;
     BuildUnipathAdjacencyHyperKmerPath( K, U, unipaths, h );
     {    double total_err = 0;
          for ( int x = 0; x < h.N( ); x++ )
          {    double in = 0, out = 0;
               for ( int j = 0; j < h.To(x).isize( ); j++ )
                    in += CN_float[ h.EdgeObjectIndexByIndexTo( x, j ) ];
               for ( int j = 0; j < h.From(x).isize( ); j++ )
                    out += CN_float[ h.EdgeObjectIndexByIndexFrom( x, j ) ];
               if ( in + out > 0 ) total_err += Abs( out - in ) / ( in + out );    }
          double conservation_err = total_err / double( h.N( ) );
          cout << "mean conservation error = " << conservation_err << endl;    }
     {    for ( int u = 0; u < nuni; u++ )
          {    int nkmers = paths[u].KmerCount( );
               if ( 10 <= nkmers && nkmers < 80 ) CN_float[u] *= 1.2;    }
          double total_err = 0;
          for ( int x = 0; x < h.N( ); x++ )
          {    double in = 0, out = 0;
               for ( int j = 0; j < h.To(x).isize( ); j++ )
                    in += CN_float[ h.EdgeObjectIndexByIndexTo( x, j ) ];
               for ( int j = 0; j < h.From(x).isize( ); j++ )
                    out += CN_float[ h.EdgeObjectIndexByIndexFrom( x, j ) ];
               if ( in + out > 0 ) total_err += Abs( out - in ) / ( in + out );    }
          double conservation_err = total_err / double( h.N( ) );
          cout << "mean conservation error = " << conservation_err << endl;    }
     return 0;
     */

     // Traverse the binary branches in the unipath graph.  Locally we require that
     // they look like this, with only two edges going out and one coming in at b:
     //
     //               ---x0--> c0
     //               | 
     //     a ---u--> b
     //               | 
     //               ---x1--> c1

     // (check to see if all vertices are distinct?)

     vec<int> to_delete;
     vec<int> bad_kills;
     for ( int u = 0; u < nuni; u++ )
     {    if ( predicted_CNs[u] > 2 ) continue;
          if ( U.From(u).size( ) != 2 ) continue;
          vec<int> x(2);
          for ( int j = 0; j < 2; j++ )
               x[j] = U.From(u)[j];
          if ( !U.To( x[0] ).solo( ) ) continue;
          if ( x[0] == u || x[1] == u ) continue;

          // For each outgoing edge xi, define its "downstream graph".  To do this
          // we successively add outgoing edges, stopping however at any vertex
          // having an incoming edge that we have not yet seen, allowing for 
          // pausing.

          vec< vec<int> > down(2);
          for ( int j = 0; j < 2; j++ )
          {    down[j].push_back( x[j] );
               vec<int> to_process, paused;
               set<int> processed;
               to_process.push_back( x[j] );
               while(1)
               {    int ndown = down[j].size( );
                    while( to_process.nonempty( ) )
                    {    int e = to_process.back( );
                         to_process.pop_back( );
                         processed.insert(e);
                         if ( U.From(e).empty( ) ) continue;
                         const vec<int>& b = U.To( U.From(e)[0] );
                         Bool outside = False;
                         for ( int z = 0; z < b.isize( ); z++ )
                         {    if ( !Member( down[j], b[z] ) )
                              {    outside = True;
                                   break;    }    }
                         if (outside) paused.push_back(e);
                         else
                         {    for ( int z = 0; z < U.From(e).isize( ); z++ )
                              {    int f = U.From(e)[z];
                                   if ( Member( processed, f ) ) continue;
                                   if ( Member( to_process, f ) ) continue;
                                   to_process.push_back(f);    
                                   down[j].push_back(f);    }    }    }
                    if ( down[j].isize( ) == ndown ) break;
                    to_process = paused;
                    paused.clear( );    }    }

          // Announce what's happening.

          if ( VERBOSITY >= 1 )
          {    cout << "\nbinary branch, edge " << u << " goes to " << x[0] 
                    << " and " << x[1] << endl;    }
          for ( int j = 0; j < 2; j++ )
          {    Sort( down[j] );
               if ( VERBOSITY >= 1 )
               {    cout << "downstream graph " << j << ":";
                    for ( int z = 0; z < down[j].isize( ); z++ )
                         cout << " " << down[j][z];
                    cout << endl;    }    }

          // Sanity check.  This used to be an assert, but evidently insanity is
          // sometimes possible.  It would be worth investigating why.

          if ( Meet( down[0], down[1] ) )
          {    if ( VERBOSITY >= 1 ) 
                    cout << "downstream graphs meet, giving up" << endl;
               continue;    }
          /* ForceAssert( !Meet( down[0], down[1] ) ); */

          // For each downstream graph, determine the minimum distance in kmers to a
          // terminal vertex.

          vec<digraph> D(2);
          vec< vec< vec<int> > > dpaths(2);
          vec<int> min_dist( 2, 1000000000 );
          for ( int j = 0; j < 2; j++ )
          {    D[j].Initialize( U, down[j] );
               D[j].AllPaths( -1, -1, dpaths[j] );
               for ( int i = 0; i < dpaths[j].isize( ); i++ )
               {    int d = 0;
                    for ( int r = 0; r < dpaths[j][i].isize( ); r++ )
                    {    int u = down[j][ dpaths[j][i][r] ];
                         d += unibases[u].isize( ) - K + 1;    }
                    min_dist[j] = Min( d, min_dist[j] );    }
               if ( VERBOSITY >= 1 )
                    cout << j << ": min kmers = " << min_dist[j] << endl;    }

          // Find eligible kmers in the downstream graphs.  

          vec< vec< pair<int,int> > > eligibles(2);
          for ( int j = 0; j < 2; j++ )
          {    vec< pair<int,int> > E;
               for ( int i = 0; i < dpaths[j].isize( ); i++ )
               {    int d = 0;
                    for ( int r = 0; r < dpaths[j][i].isize( ); r++ )
                    {    int u = down[j][ dpaths[j][i][r] ];
                         int elig 
                              = Min( unibases[u].isize( ) - K + 1, min_dist[j] - d );
                         E.push( u, elig );
                         d += unibases[u].isize( ) - K + 1;
                         if ( d >= min_dist[j] ) break;    }    }    
               UniqueSort(E);
               for ( int r = 0; r < E.isize( ); r++ )
               {    int u = E[r].first, s;
                    for ( s = r + 1; s < E.isize( ); s++ )
                         if ( E[s].first != u ) break;
                    if ( VERBOSITY >= 1 )
                    {    cout << j << ": eligible " << u << ".1-" 
                              << E[s-1].second << endl;    }
                    eligibles[j].push( u, E[s-1].second );
                    r = s - 1;    }    }

          // For each eligible kmer, define its start position relative to
          // the start of the branch.  We track this by kmer numbers.

          vec< vec<longlong> > kmer_number(2);
          vec< vec<int> > kmer_number_start(2);
          for ( int j = 0; j < 2; j++ )
          {    for ( int r = 0; r < eligibles[j].isize( ); r++ )
               {    int w = eligibles[j][r].first;

                    // Track back to x[j] easily, or give up.

                    int w_start_rel_xj = 0, y = w;
                    while(1)
                    {    if ( y == x[j] || U.To(y).isize( ) != 1 ) break;
                         y = U.To(y)[0];
                         w_start_rel_xj += unibases[y].isize( ) - K + 1;    }
                    if ( y != x[j] ) continue;

                    // Go through the kmer numbers.

                    const KmerPath& p = paths[w];
                    int n = eligibles[j][r].second;
                    KmerPathLoc loc = p.Begin( );
                    for ( int s = 0; s < n; s++ )
                    {    kmer_number[j].push_back( loc.GetKmer( ) );
                         kmer_number_start[j].push_back( w_start_rel_xj + s );
                         loc += 1;    }    }
               SortSync( kmer_number[j], kmer_number_start[j] );    }

          // Get reads for each branch.  We also track their orientations.

          vec< vec<int> > ids(2);
          vec< vec<Bool> > fw(2);
          for ( int j = 0; j < 2; j++ )
          {    for ( int r = 0; r < eligibles[j].isize( ); r++ )
               {    int u = eligibles[j][r].first, n = eligibles[j][r].second;
                    int pos = 0;
                    for ( int s = 0; s < paths[u].NSegments( ); s++ )
                    {    longlong start = paths[u].Start(s), stop = paths[u].Stop(s);
                         if ( pos + stop - start + 1 > n )
                              stop -= ( pos + stop - start + 1 ) - n;
                         pos += stop - start + 1;
                         vec<longlong> con;
                         Contains( pathsdb, KmerPathInterval( start, stop ), con );
                         for ( int l = 0; l < con.isize( ); l++ )
                         {    const tagged_rpint& t = pathsdb[ con[l] ];
                              int id = t.ReadId( );
                              if ( id >= (int) nuni ) 
                              {    ids[j].push_back( id-nuni );
                                   fw[j].push_back( t.Fw( ) );    }    }
                         if ( pos == n ) break;    }    }

               // Note that the following effectively picks an orientation at
               // random if a read is assigned two orientations.  This is not
               // correct.

               UniqueSortSync( ids[j], fw[j] );    }

          // Remove reads that are shared by the branches.

          vec<int> ids_share = Intersection( ids[0], ids[1] );
          for ( int j = 0; j < 2; j++ )
          {    vec<Bool> to_delete( ids[j].size( ), False );
               for ( int i = 0; i < ids[j].isize( ); i++ )
                    if ( BinMember( ids_share, ids[j][i] ) ) to_delete[i] = True;
               EraseIf( ids[j], to_delete ), EraseIf( fw[j], to_delete );
               cout << j << ": " << ids[j].size( ) << " reads hit it" << endl;
               if (PRINT_READ_IDS)
               {    for ( int l = 0; l < ids[j].isize( ); l++ )
                    {    if ( l > 0 ) cout << ",";
                         cout << ids[j][l];    }
                    cout << "\n";    }    }

          // Remove duplicate pairs.  Note that this would not be needed if 
          // duplicate fragment pairs were removed early in the process.  Duplicates
          // are defined by agreement along the first 20 bases.

          const int start_bases = 20;
          for ( int j = 0; j < 2; j++ )
          {    vec<Bool> dup( ids[j].size( ), False );
               vec< triple<basevector,basevector,int> > starts;
               for ( int i = 0; i < ids[j].isize( ); i++ )
               {    int64_t id1 = ids[j][i];
                    if ( id1 >= nfrags ) continue;
                    int64_t id2 = pairs.getPartnerID(id1);
                    if ( id2 >= 0 )
                    {    const basevector &r1 = all[nuni+id1], &r2 = all[nuni+id2];
                         if ( r1.isize( ) >= start_bases 
                              && r2.isize( ) >= start_bases )
                         {    basevector s1(r1, 0, start_bases); 
                              basevector s2(r2, 0, start_bases);
                              starts.push( s1, s2, id1 );    }    }    }
               Sort(starts);
               for ( int l1 = 0; l1 < starts.isize( ); l1++ )
               {    int l2;
                    for ( l2 = l1 + 1; l2 < starts.isize( ); l2++ )
                    {    if ( starts[l2].first != starts[l1].first ) break;
                         if ( starts[l2].second != starts[l1].second ) break;    }
                    for ( int m = l1 + 1; m < l2; m++ )
                    {    int p = BinPosition( ids[j], starts[m].third );
                         if ( p >= 0 ) dup[p] = True;    }
                    l1 = l2 - 1;    }
               cout << j << ": deleting " << Sum(dup) << " duplicate reads" << "\n";
               EraseIf( ids[j], dup ), EraseIf( fw[j], dup );    }

          // For each read, determine its start point relative to the origin of the
          // branch.  For the moment we do not worry about multiple start points for
          // a given read.  Note that we don't check to see that start positions for
          // all reads are assigned.

          vec< vec<int> > start(2);
          for ( int j = 0; j < 2; j++ )
          {    start[j].resize( ids[j].size( ), 0 );
               for ( int i = 0; i < ids[j].isize( ); i++ )
               {    int id = ids[j][i];
                    const KmerPath& p = ( fw[j][i] ? paths : pathsrc )[nuni+id];
                    KmerPathLoc l = p.Begin( );
                    int n = p.KmerCount( );
                    for ( int r = 0; r < n; r++ )
                    {    longlong x = l.GetKmer( );
                         int low = lower_bound( kmer_number[j].begin( ),
                              kmer_number[j].end( ), x ) - kmer_number[j].begin( );
                         int high = upper_bound( kmer_number[j].begin( ),
                              kmer_number[j].end( ), x ) - kmer_number[j].begin( );
                         for ( int z = low; z < high; z++ )
                         {    int s = kmer_number_start[j][z] - r;
                              start[j][i] = s; // note could assign different values!
                                   }
                         l += 1;    }    }    }

          // Define counts.

          int x1 = x[0], x2 = x[1];
          int md = Min( min_dist[0], min_dist[1] );
          double mult1 = double(md) / double( min_dist[0] );
          double mult2 = double(md) / double( min_dist[1] );
          long double h1 = (long double)
               ( ids[0].size( ) + extra_hits[x1] + extra[x1] ) * mult1;
          long double h2 = (long double)
               ( ids[1].size( ) + extra_hits[x2] + extra[x2] ) * mult2;

          // Compare x1 to x2.  We define four categories:
          // (1) indel of up to 4 bases;
          // (2) mismatch adjacent to 10 or more copies of a motif of size up to 4;
          // (3) some other mismatch;
          // (4) something else.

          const int min_motif_copies = 10;
          const int max_motif_size = 4;
          Bool simple = False, mismatch = False, simple_repeat_mismatch = False;
          const basevector &X1 = unibases[x1], &X2 = unibases[x2];
          int n1 = unibases[x1].size( ), n2 = unibases[x2].size( );
          Bool agree = True;
          for ( int j = K; j < Min( n1, n2 ); j++ )
          {    if ( X1[j] != X2[j] ) 
               {    agree = False;
                    break;    }    }
          if (agree) 
          {    simple = True;
               cout << "x1 and x2 differ by a mismatch\n";
               mismatch = True;
               for ( int pass = 1; pass <= 2; pass++ )
               {    for ( int m = 1; m <= max_motif_size; m++ )
                    {    if ( pass == 1 && (K-1) - min_motif_copies*m < 0 )
                              continue;
                         if ( pass == 2 
                              && (K-1) + min_motif_copies*m >= Min( n1, n2 ) )
                         {    continue;    }
                         Bool fail = False;
                         for ( int l = 0; l < min_motif_copies - 1; l++ )
                         {    for ( int r = 0; r < m; r++ )
                              {    int p;
                                   if ( pass == 1 ) // to left of mismatch
                                   {    p = (K-1) - (m*l+r+1) - m;    }
                                   else // to right of mismatch
                                   {    p = (K-1) + (m*l+r+1);    }
                                   if ( X1[p] != X1[p+m] )
                                   {    fail = True;
                                        break;    }    }
                              if (fail) break;    }
                         if ( !fail ) 
                         {    simple_repeat_mismatch = True;
                              cout << "x1 and x2 differ by a simple "
                                   << "repeat mismatch\n";
                              break;    }    }
                    if (simple_repeat_mismatch) break;    }    }
          else
          {    int in;
               for ( in = -max_motif_size; in <= max_motif_size; in++ )
               {    Bool agree = True;
                    if ( in > 0 )
                    {    for ( int j = K; j < Min( n1-in, n2 ); j++ )
                         {    if ( X1[j+in] != X2[j] ) 
                              {    agree = False; break;    }    }
                         if (agree) 
                         {    simple = True; break;    }    }
                    if ( in < 0 )
                    {    for ( int j = K; j < Min( n1, n2-in ); j++ )
                         {    if ( X1[j] != X2[j+in] ) 
                              {    agree = False; break;    }    }
                         if (agree) 
                         {    simple = True; break;    }    }    }
               if (simple)
               {    cout << "x1 and x2 differ by an indel of size " << Abs(in)
                         << "\n";    }    }

          // In case (4), i.e. neither a mismatch nor an indel, if there are two 
          // reads supporting each branch, we give up.

          if ( h1 > 1 && h2 > 1 && !simple )
          {    cout << "Not a simple difference between x1 and x2, and multi-read "
                    << "support,\n" << "so not considering either for deletion.\n";
               continue;    }

          // Look at quality scores in the mismatch case.

          vec< vec<int> > qs(2);
          if ( /* mismatch */ !simple_repeat_mismatch )
          {    int raw1 = ids[0].size( ) + extra_hits[x1];
               int raw2 = ids[1].size( ) + extra_hits[x2];
               char b1 = unibases[x1][K-1], b2 = unibases[x2][K-1];
               for ( int j = 0; j < 2; j++ )
               {    for ( int l = 0; l < ids[j].isize( ); l++ )
                    {    
                         /*
                         basevector b = all[ nuni + ids[j][l] ];
                         qualvector q = quals[ ids[j][l] ];
                         */

                         basevector b;
                         qualvector q;
                         if ( ids[j][l] < nfrags )
                         {    b = bases0[ ids[j][l] ];
                              q = quals0[ ids[j][l] ];    }
                         else
                         {    b = jbases0[ ids[j][l] - nfrags ];
                              q = jquals0[ ids[j][l] - nfrags ];    }

                         if ( !fw[j][l] ) 
                         {    b.ReverseComplement( );
                              q.ReverseMe( );    }
                         int p = K - 1 - start[j][l];
                         if ( p < 0 || p >= b.isize( ) )
                         {    if ( VERBOSITY >= 2 )
                              {    cout << j << ": read " << ids[j][l]
                                        << " does not report\n";    }    }
                         else
                         {    if ( VERBOSITY >= 1 )
                              {    cout << j << ": read " << ids[j][l] 
                                        << " reports " << as_base( b[p] ) << "[" 
                                        << int(q[p]) << "]\n";    }
                              if ( j == 0 && b[p] == b2 )
                              {    --raw1;
                                   ++raw2;
                                   qs[1].push_back( q[p] );    }
                              else if ( j == 1 && b[p] == b1 )
                              {    --raw2;
                                   ++raw1;
                                   qs[0].push_back( q[p] );    }
                              else qs[j].push_back( q[p] );    }    }
                    ReverseSort( qs[j] );    }    
               h1 = (long double)( raw1 ) * mult1;
               h2 = (long double)( raw2 ) * mult2;    }

          // Pick branch.

          if ( h1 <= 1 && h2 <= 1 ) continue;
          Bool swapped = False;
          if ( h2 < h1 ) 
          {    swap( h1, h2 );
               swap( x1, x2 );    
               swap( qs[0], qs[1] );
               swapped = True;    }
          if ( h1 < 1 ) h1 = 1;
          long double n = h1 + h2;

          // In the case of a mismatch that is not a simple repeat mismatch, if 
          // there are two reads at Q30+ supporting the base, give up.

          if ( /* mismatch && */ !simple_repeat_mismatch )
          {    if ( qs[0].size( ) >= 2 && qs[0][1] >= 30 )
               {    cout << "Mismatch, but not a simple repeat mismatch between "
                         << "x1 and x2,\n" << "and have double Q30 support, "
                         << "so not considering deletion.\n";
                    continue;    }
               if ( Sum( qs[0] ) >= 100 )
               {    cout << "Mismatch, but not a simple repeat mismatch between "
                         << "x1 and x2,\n" << "and have Q support of at least 100, "
                         << "so not considering deletion.\n";
                    continue;    }
               if ( Sum( qs[0] ) >= Sum( qs[1] ) )
               {    cout << "Mismatch, but not a simple repeat mismatch between "
                         << "x1 and x2,\n" << "and Q support inadequate to "
                         << "support deletion.\n";
                    continue;    }    }

          // Assume that p(error) = h1/n, in the haploid case, or 0.5 in the diploid
          // case.  We require that the probability of h2 or more errors (out of n) 
          // is <= 0.001.  Note that the output of BinomialSum might be wrong if its
          // first parameter is too large, but it is still hopefully the case that 
          // the test gives the right answer.

          long double p = ( predicted_CNs[u] <= 1 ? h1/n : 0.5 );
          long double p_wild 
               = 1.0 - BinomialSum( int(round(n)), int(round(h2-1)), p );
          if ( p_wild > max_p_opp ) continue;
          if ( VERBOSITY >= 1 )
          {    cout << "u = " << u << ", x1 = " << ( !swapped ? x1 : x2 )
                    << ", x2 = " << ( !swapped ? x2 : x1 ) << ", h1 = " 
                    << ( !swapped ? h1 : h2 ) << ", h2 = " << ( !swapped ? h2 : h1 )
                    << ", p_wild = " << p_wild;    }
          if ( true_CN[x1] > 0 ) 
          {    if ( VERBOSITY >= 1 ) cout << " ERROR!";
               bad_kills.push_back( x1, to_rc[x1] );    }
          if ( VERBOSITY >= 1 )
          {    cout << "\n";
               if (mismatch) 
               {    if ( !swapped ) 
                         cout << "q1 = " << Sum(qs[0]) << ", q2 = " << Sum(qs[1]);
                    else cout << "q1 = " << Sum(qs[1]) << ", q2 = " << Sum(qs[0]);
                    cout << "\n";    }
               cout << "recommend deleting unipath " << x1 
                    << " in favor of " << x2 << "\n";    
               cout << "recommend deleting unipath " << to_rc[x1] << "\n";    }
          to_delete.push_back( x1, to_rc[x1] );    }

     // Report results.

     UniqueSort(to_delete), UniqueSort(bad_kills);
     if ( !VALIDATE )
          cout << "\n" << to_delete.size( ) << " unipaths deleted" << endl;
     else
     {    cout << "\n" << to_delete.size( ) - bad_kills.size( ) 
               << " bad unipaths deleted" << endl;
          cout << bad_kills.size( ) << " good unipaths deleted" << endl;    }

     // Regenerate unipaths.  This presumably not the right solution.

     {
     cout << Date( ) << ": " << unibases.size( ) << " old unibases" << endl;
     vec<Bool> uni_to_delete( nuni, False );
     for ( int i = 0; i < to_delete.isize( ); i++ )
          uni_to_delete[ to_delete[i] ] = True;
     unibases.EraseIf(uni_to_delete);
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     size_t nuni = unibases.size( );
     for ( size_t u = 0; u < nuni; u++ )
     {    for ( int j = 0; j < nexts[u].isize( ); j++ )
          {    basevector b1 = unibases[u], b2 = unibases[ nexts[u][j] ];
               b1.SetToSubOf( b1, b1.isize( ) - K, K );
               b2.SetToSubOf( b2, K-1, 1 );
               basevector b = Cat( b1, b2 );
               unibases.push_back_reserve(b);    }    }
     size_t N = unibases.size( );
     vecKmerPath paths, pathsrc, unipaths;
     vec<tagged_rpint> pathsdb, unipathsdb;
     const int threads_to_use = 2 * NUM_THREADS / 3;
     Mkdir777( run_dir + "/tmp" );
     ReadsToPathsCoreY( unibases, K, paths, pathsrc, 
          pathsdb, run_dir + "/tmp", threads_to_use );
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
     cout << Date( ) << ": " << unibases.size( ) << " new unibases" << endl;

     // Compute summary stats from unibases.

     if (STATS)
     {    cout << Date( ) << ": computing summary stats from new unibases" 
               << endl;
          vecbasevector genome( data_dir + "/genome.fastb" );
          vec<basevector> missing_genomic_kmers;
          int64_t non_genomic_kmers = 0, N50_unibase = 0;
          if ( K == 96 )
          {    UnibaseSummaryStats<96>( unibases, genome, missing_genomic_kmers, 
                    non_genomic_kmers, N50_unibase );    }
          else ForceAssert( 0 == 1 );
          cout << "missing genomic kmers = " << missing_genomic_kmers.size( ) 
               << endl;
          cout << "non-genomic kmers = " << non_genomic_kmers << endl;
          cout << "N50 unibase = " << N50_unibase << endl;    }
     
     // Write new files.

     if (WRITE)
     {    String KS = ToString(K);
          unibases.WriteAll( run_dir + "/" + UNIBASES_OUT + ".unibases.k" + KS );
          paths.WriteAll( run_dir + "/" + UNIBASES_OUT + ".paths.k" + KS );
          pathsrc.WriteAll( run_dir + "/" + UNIBASES_OUT + ".paths_rc.k" + KS );
          unipaths.WriteAll( run_dir + "/" + UNIBASES_OUT + ".unipaths.k" + KS );
          BinaryWrite3( run_dir + "/" + UNIBASES_OUT + ".pathsdb.k" + KS, pathsdb );
          BinaryWrite3( run_dir + "/" + UNIBASES_OUT + ".unipathsdb.k" + KS, 
               unipathsdb );    }    }    }
