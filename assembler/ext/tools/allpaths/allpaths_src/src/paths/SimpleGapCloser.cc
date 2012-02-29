///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2012) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SimpleGapCloser.  Look for data
//
//                          ----------read-----------
//       -------u1--------------                --------------u2-----------
//
// in which unibases u1 [resp. u2] are terminal on the right [resp. left], and
// the read perfectly matches both the last 20 bases of u1 and the first 20 bases
// of u2.  These reads are treated as possible patches for a putative gap between
// u1 and u2.  Other possible placements for the read are also considered, and the
// patch is only allowed if the given placement is best.
//
// The above explanation was written before several changes were made to the
// algorithm.  It is not up to date.
//
// An important point is that the algorithm does not use pairing information, and
// thus may work for 'crappy' libraries when other unipath gap closing methods fail.

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS
// MakeDepend: dependency QueryLookupTable

#include <omp.h>

#include "Basevector.h"
#include "MainTools.h"
#include "ParallelVecUtilities.h"
#include "Qualvector.h"
#include "kmers/KmerRecord.h"
#include "lookup/LookAlign.h"
#include "paths/GetNexts.h"
#include "paths/KmerBaseBroker.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/RemodelGapTools.h"
#include "paths/UnibaseUtils.h"
#include "paths/Unipath.h"

class bridge {

     public:

     bridge( ) { }
     bridge( const int uid1, const int uid2, const int gap, 
          const basevector& bases ) : uid1(uid1), uid2(uid2), gap(gap), 
          bases(bases) { }

     friend Bool operator<( const bridge& b1, const bridge& b2 )
     {    if ( b1.uid1 < b2.uid1 ) return True;
          if ( b1.uid1 > b2.uid1 ) return False;
          if ( b2.uid2 < b2.uid2 ) return True;
          if ( b2.uid2 > b2.uid2 ) return False;
          if ( b1.gap < b2.gap ) return True;
          if ( b1.gap > b2.gap ) return False;
          if ( b1.bases < b2.bases ) return True;
          return False;    }

     int uid1, uid2;
     int gap;
     basevector bases;

};

void QLT( const basevector& b1, const basevector& b2, const basevector& b3,
     const String& data_dir, const String& run_dir )
{    vecbasevector B;
     B.push_back(b1), B.push_back(b2), B.push_back(b3);
     String fn = run_dir + "/SimpleGapCloser.QLT.tmp.fastb";
     B.WriteAll(fn);
     vec<String> out = AllOfOutput( "QueryLookupTable K=12 MM=12 MC=0.15 "
          "SEQS=" + fn + " L=" + data_dir + "/genome.lookup "
          + "VISUAL=True NH=True QUIET=True PARSEABLE=True" );
     for ( int i = 0; i < out.isize( ); i++ )
     {    if ( out[i].Contains( "2fw", 0 ) || out[i].Contains( "2rc", 0 ) )
          {    if ( out[i].Contains( " 0 mismatches/0 indels" ) )
               {    int len = out[i].Between( "of ", ")" ).Int( );
                    if ( out[i].Contains( "from 0-" + ToString(len) + " " ) )
                    {    cout << "\nperfect\n";
                         return;    }    }    }    }
     for ( int i = 0; i < out.isize( ); i++ )
          if ( !out[i].Contains( "QUERY", 0 ) ) cout << out[i] << "\n";    }

int main(int argc, char *argv[])
{
     RunTime( );

     BeginCommandArguments;
     CommandArgument_String(PRE);
     CommandArgument_String(DATA);
     CommandArgument_String(RUN);
     CommandArgument_Int(K);
     CommandArgument_String(HEAD_IN);
     CommandArgument_String(HEAD_OUT);
     CommandArgument_Bool_OrDefault(WRITE, True);
     CommandArgument_Bool_OrDefault(VALIDATE, False);
     CommandArgument_Bool_OrDefault(VERBOSE, False);
     CommandArgument_Int_OrDefault_Doc(NUM_THREADS, -1,
          "number of threads to be used (use all available if negative)");
     EndCommandArguments;

     // Thread control

     NUM_THREADS = configNumThreads(NUM_THREADS);
     omp_set_num_threads( NUM_THREADS );

     // Heuristics.

     const int max_placements = 10;
     const int max_paths = 100;
     const int max_pushes = 1000;
     const int top_qtotal = 50;
     const int top_maxq = 20;
     const int min_bridges = 2;
     const int min_perf_ext = 15;

     // Define directories, etc.

     String data_dir = PRE + "/" + DATA; 
     String run_dir = PRE + "/" + DATA + "/" + RUN;

     // Load unibases and find terminal ones.

     cout << Date( ) << ": loading unibases" << endl;
     String unibases_file = run_dir + "/" + HEAD_IN + ".unibases.k" + ToString(K);
     vecbasevector unibases(unibases_file);
     vec< vec<int> > nexts;
     GetNexts( K, unibases, nexts );
     vec<int> to_rc;
     UnibaseInvolution( unibases, to_rc, K );
     vec<int> starts, stops;
     for ( size_t u = 0; u < unibases.size( ); u++ )
     {    if ( nexts[u].empty( ) ) starts.push_back(u);
          if ( nexts[ to_rc[u] ].empty( ) ) stops.push_back(u);    }

     // Define the terminal kmers of the starts and stops.

     const int L = 20;
     vec< kmer<L> > startsb( starts.size( ) ), stopsb( stops.size( ) );
     for ( int i = 0; i < starts.isize( ); i++ )
     {    const basevector& U = unibases[ starts[i] ];
          startsb[i].SetToSubOf( U, U.isize( ) - L );    }
     for ( int i = 0; i < stops.isize( ); i++ )
     {    const basevector& U = unibases[ stops[i] ];
          stopsb[i].SetToSubOf( U, 0 );    }
     vec<int> startsi(starts), stopsi(stops);
     SortSync( startsb, startsi ), SortSync( stopsb, stopsi );

     // Create lookup table for kmers.

     vec< triple<kmer<L>,int,int> > kmers_plus;
     MakeKmerLookup0( unibases, kmers_plus );
     vec< kmer<L> > kmers( kmers_plus.size( ) );
     for ( size_t i = 0; i < kmers.size( ); i++ )
          kmers[i] = kmers_plus[i].first;

     // Scan the reads to find ones that contain start and stop kmers.

     vec<bridge> bridges;
     for ( int pass = 1; pass <= 2; pass++ )
     {    String lib = ( pass == 1 ? "jump" : "frag" );
          cout << Date( ) << ": loading " << lib << " reads" << endl;
          vecbasevector bases( run_dir + "/" + lib + "_reads_filt.fastb" );
          vecqualvector quals( run_dir + "/" + lib + "_reads_filt.qualb" );
          cout << Date( ) << ": finding extensions" << endl;
          const int batch = 10000;
          #pragma omp parallel for
          for ( size_t idx = 0; idx < bases.size( ); idx += batch )
          {    kmer<L> x;
               vec<bridge> bridges_this;
               for ( size_t id = idx; id < Min( bases.size( ), idx + batch ); id++ )
               {    const basevector& b = bases[id];
     
                    // First see if the read could bridge a 'gap' between two 
                    // unibases.

                    vec< pair<int,int> > ext_right, ext_left;
                    for ( int j = 0; j <= b.isize( ) - L; j++ )
                    {    x.SetToSubOf( b, j );
                         int low = LowerBound(startsb, x); 
                         int high = UpperBound(startsb, x);
                         for ( int i = low; i < high; i++ )
                              ext_right.push( j, startsi[i] );
                         low = LowerBound(stopsb, x), high = UpperBound(stopsb, x);
                         for ( int i = low; i < high; i++ )
                              ext_left.push( j, stopsi[i] );    }    
                    if ( ext_left.empty( ) || ext_right.empty( ) ) continue;

                    // Define bridges.  Test for implied overlap.  We require that
                    // there is exactly one bridge.

                    vec<bridge> bridges_this_this;
                    vec< pair<int,int> > vx;
                    Bool too_many = False;
                    int leftx = 0, rightx = 0;
                    for ( int j1 = 0; j1 < ext_right.isize( ); j1++ )
                    {    if (too_many) break;
                         for ( int j2 = 0; j2 < ext_left.isize( ); j2++ )
                         {    int u1 = ext_right[j1].second; 
                              int u2 = ext_left[j2].second;
                              int start = ext_right[j1].first + L;
                              int stop = ext_left[j2].first;
                              const basevector &U1 = unibases[u1]; 
                              const basevector &U2 = unibases[u2];
                              int overlap = start - stop;
                              Bool OK = True;
                              if ( overlap > 0 )
                              {    for ( int j = stop; j < start; j++ )
                                   {    if ( U1[ U1.isize( ) - (j-stop+1) ]
                                             != U2[ overlap - (j-stop+1) ] )
                                        {    OK = False;
                                             break;    }    }    }
                              if ( !OK || overlap == U1.isize( ) ) continue;
                              if ( bridges_this_this.nonempty( ) )
                              {    too_many = True;
                                   break;    }
                              leftx = ext_right[j1].first;
                              rightx = ext_left[j2].first;
                              vx.push( u1, ext_right[j1].first 
                                   - ( unibases[u1].isize( ) - L ) );
                              vx.push( u2, ext_left[j2].first );
                              basevector b;
                              if ( overlap < 0 )
                                   b.SetToSubOf( bases[id], start, stop-start );
                              bridges_this_this.push( u1, u2, stop-start, b );    
                              b.ReverseComplement( );
                              bridges_this_this.push( to_rc[u2], to_rc[u1], 
                                   stop-start, b );    }    }    
                    if ( too_many || bridges_this_this.size( ) != 2 ) continue;

                    // Check for the required perfect extension.

                    int start = leftx + L, stop = rightx;
                    int overlap = start - stop;
                    if ( overlap > 0 )
                    {    int u1 = vx[0].first, u2 = vx[1].first;
                         int perf1 = 0, perf2 = 0, max_perf1 = 0, max_perf2 = 0;
                         for ( int rpos = 0; rpos < vx[1].second; rpos++ )
                         {    int u1pos = rpos - vx[0].second;
                              if ( u1pos < 0 || u1pos >= unibases[u1].isize( ) ) 
                                   continue;
                              if ( bases[id][rpos] == unibases[u1][u1pos] ) perf1++;
                              else 
                              {    max_perf1 = Max( max_perf1, perf1 );
                                   perf1 = 0;    }    }
                         max_perf1 = Max( max_perf1, perf1 );
                         for ( int rpos = vx[0].second + unibases[u1].isize( );
                         rpos < bases[id].isize( ); rpos++ )
                         {    int u2pos = rpos - vx[1].second;
                              if ( u2pos < 0 || u2pos >= unibases[u2].isize( ) ) 
                                   continue;
                              if ( bases[id][rpos] == unibases[u2][u2pos] ) perf2++;
                              else 
                              {    max_perf2 = Max( max_perf2, perf2 );
                                   perf2 = 0;    }    }
                         max_perf2 = Max( max_perf2, perf2 );
                         if ( max_perf1 < min_perf_ext || max_perf2 < min_perf_ext )
                              continue;    }

                    // Give up if the read is subsumed by either unipath.

                    /*
                    Bool subsumed = False;
                    for ( int j = 0; j < 2; j++ )
                    {    if ( vx[j].second <= 0 && bases[id].isize( ) 
                              <= vx[j].second + unibases[ vx[j].first ].isize( ) )
                         {    subsumed = True;    }    }
                    if (subsumed) continue;
                    */

                    // Now find all placements for the read.

                    vec< pair<int,int> > v;
                    for ( int j = 0; j <= b.isize( ) - L; j++ )
                    {    x.SetToSubOf( b, j );
                         int64_t low = LowerBound(kmers, x); 
                         int64_t high = UpperBound(kmers, x);
                         for ( int64_t l = low; l < high; l++ )
                         {    int u = kmers_plus[l].second; 
                              int upos = kmers_plus[l].third;
                              v.push( u, j - upos );    }    }
                    UniqueSort(v);
                    int p1 = BinPosition( v, vx[0] ), p2 = BinPosition( v, vx[1] );
                    ForceAssert( p1 >= 0 ); ForceAssert( p2 >= 0 );
                    if ( v.isize( ) > max_placements ) continue;

                    // Form the placements into a digraph.

                    int N = v.size( );
                    vec< vec<int> > from(N), to(N);
                    for ( int j1 = 0; j1 < N; j1++ )
                    for ( int j2 = 0; j2 < N; j2++ )
                    {    int u1 = v[j1].first, u2 = v[j2].first;
                         const basevector &U1 = unibases[u1], &U2 = unibases[u2];
                         int o1 = v[j1].second, o2 = v[j2].second;
                         if ( !( o1 < o2 ) ) continue;
                         int stop1 = o1 + U1.isize( ), start2 = o2;
                         int overlap = stop1 - start2;
                         Bool mismatch = False;
                         if ( overlap > U2.isize( ) ) continue;
                         for ( int j = 0; j < overlap; j++ )
                         {    if ( U1[ U1.isize( ) - j - 1 ]
                                   !=  U2[ overlap - j - 1 ] )
                              {    mismatch = True;
                                   break;    }    }
                         if (mismatch) continue;
                         from[j1].push_back(j2), to[j2].push_back(j1);    }
                    for ( int j = 0; j < N; j++ )
                    {    Sort( from[j] ), Sort( to[j] );    }
                    digraph G( from, to );
                    vec< vec<int> > paths;
                    if ( !G.AllPaths( -1, -1, paths, max_paths, False, max_pushes ) )
                         continue;

                    // Put the 'standard path' in front.

                    vec<int> pcore;
                    pcore.push_back( p1, p2 );
                    int pi = Position( paths, pcore );
                    if ( pi >= 0 ) 
                    {    if ( pi != 0 ) swap( paths[0], paths[pi] );    }
                    else paths.push_front(pcore);

                    // Check its quality scores.

                    const vec<int>& p = paths[0];
                    vec<int> errs;
                    for ( int l = 0; l < p.isize( ); l++ )
                    {    int u = v[ p[l] ].first, pos = v[ p[l] ].second;
                         const basevector &U = unibases[u], &R = bases[id];
                         int ulen = U.size( ), rlen = R.size( );
                         int start = Max(0, pos), stop = Min(rlen, pos + ulen);
                         for ( int upos = 0; upos < U.isize( ); upos++ )
                         {    int rpos = upos + pos;
                              if ( rpos < 0 ) continue;
                              if ( rpos >= rlen ) break;
                              if ( R[rpos] != U[upos] ) 
                                   errs.push_back(rpos);    }    }
                    UniqueSort(errs);
                    int q = 0, maxq = 0;
                    for ( int l = 0; l < errs.isize( ); l++ )
                    {    q += quals[id][ errs[l] ];
                         maxq = Max( maxq, (int) quals[id][ errs[l] ] );    }
                    if ( q > top_qtotal || maxq > top_maxq ) continue;

                    // Compute the multiplicity of the end kmers.

                    kmer<L> xleft, xright;
                    xleft.SetToSubOf( bases[id], leftx );
                    xright.SetToSubOf( bases[id], rightx );
                    int64_t nleft 
                         = UpperBound(kmers, xleft) - LowerBound(kmers, xleft); 
                    int64_t nright 
                         = UpperBound(kmers, xright) - LowerBound(kmers, xright); 

                    // Require uniqueness.

                    if ( overlap < L && ( nleft > 1 || nright > 1 ) ) continue;
                    if ( nleft > 2 || nright > 2 ) continue;

                    // Save bridges and print info.

                    bridges_this.append(bridges_this_this);
                    if ( !VERBOSE ) continue;
                    #pragma omp critical
                    {    cout << "\n=== jump read " << id << " ===" << endl;
                         PRINT2( nleft, nright );
                         for ( int j = 0; j < paths.isize( ); j++ )
                         {    const vec<int>& p = paths[j];
                              vec<int> cov, errs;
                              cout << "\nsee path:\n";
                              for ( int l = 0; l < p.isize( ); l++ )
                              {    int u = v[ p[l] ].first, pos = v[ p[l] ].second;
                                   const basevector &U = unibases[u], &R = bases[id];
                                   int ulen = U.size( ), rlen = R.size( );
                                   int start = Max(0, pos); 
                                   int stop = Min(rlen, pos + ulen);
                                   cout << u << ": " << start << " to " << stop
                                        << ( stop == rlen ? " = end" : "" ) << "\n";
                                   for (int upos = 0; upos < U.isize( ); upos++)
                                   {    int rpos = upos + pos;
                                        if ( rpos < 0 ) continue;
                                        if ( rpos >= rlen ) break;
                                        cov.push_back(rpos);
                                        if ( R[rpos] != U[upos] ) 
                                             errs.push_back(rpos);    }    }
                              UniqueSort(cov), UniqueSort(errs);
                              cout << "errs:";
                              for ( int l = 0; l < errs.isize( ); l++ )
                                   cout << " " << errs[l];
                              int q = 0, maxq = 0;
                              for ( int l = 0; l < errs.isize( ); l++ )
                              {    q += quals[id][ errs[l] ];
                                   maxq = Max( maxq, 
                                        (int) quals[id][ errs[l] ] );    }
                              int ncov = cov.size( );
                              cout << "\ncov = " << ncov << ", qual = " << q 
                                   << ", " << "qual/cov = " << double(q)/double(ncov)
                                   << ", maxq = " << maxq << "\n";    }
                         {    cout << "\nbridges" << endl;
                              for ( int j1 = 0; j1 < ext_right.isize( ); j1++ )
                              {    for (int j2 = 0; j2 < ext_left.isize(); j2++)
                                   {    int u1 = ext_right[j1].second; 
                                        int u2 = ext_left[j2].second;
                                        int start = ext_right[j1].first + L;
                                        int stop = ext_left[j2].first;
                                        PRINT4( u1, u2, start, stop );    
                                             }    }    }    }    }
               #pragma omp critical
               {    bridges.append(bridges_this);    }    }    }

     // Create patches.

     vecbasevector merged_paths;
     int bridge_count = 0;
     Sort(bridges);
     for ( int j = 0; j < bridges.isize( ); j++ )
     {    const bridge& b = bridges[j];
          int k, u1 = b.uid1, u2 = b.uid2;
          for ( k = j+1; k < bridges.isize( ); k++ )
               if ( bridges[k].uid1 != u1 || bridges[k].uid2 != u2 ) break;
          vec<int> overlaps;
          vec<basevector> gaps;
          for ( int l = j; l < k; l++ )
          {    if ( bridges[l].gap < 0 ) overlaps.push_back( -bridges[l].gap );
               else gaps.push_back( b.bases );    }
          Sort(overlaps), Sort(gaps);
          const basevector &U1 = unibases[u1], &U2 = unibases[u2];
          int n1 = U1.size( ), n2 = U2.size( );
          Bool bridged = False;
          ostringstream msg;
          msg << "\nbridges from " << u1 << "[l=" << n1 << "] to " << u2 << "[l=" 
               << n2 << "] (" << to_rc[u2] << " to " << to_rc[u1] << "):\n";
          for ( int r = 0; r < overlaps.isize( ); r++ )
          {    int s = overlaps.NextDiff(r);
               if ( s-r >= min_bridges )
               {    if (VERBOSE)
                    {    cout << ( bridged ? "" : msg.str( ) ) << "[" << s-r 
                              << "] overlap " << overlaps[r] << endl;    }
                    bridged = True;
                    basevector M = Cat( basevector(U1, 0, n1 - overlaps[r]), U2 );
                    merged_paths.push_back_reserve(M);
                    if (VALIDATE) QLT( U1, U2, M, data_dir, run_dir );    }
               r = s - 1;    }
          for ( int r = 0; r < gaps.isize( ); r++ )
          {    int s = gaps.NextDiff(r);
               if ( s-r >= min_bridges )
               {    if (VERBOSE)
                    {    cout << ( bridged ? "" : msg.str( ) ) << "[" << s-r 
                              << "] gap " << gaps[r].ToString( ) << endl;    }
                    bridged = True;
                    basevector M = Cat( U1, gaps[r], U2 );
                    merged_paths.push_back_reserve(M);
                    if (VALIDATE) QLT( U1, U2, M, data_dir, run_dir );    }
               r = s - 1;    }
          if (bridged) bridge_count++;
          j = k - 1;    }
     if (VERBOSE) cout << "\n";
     cout << Date( ) << ": bridged " << bridge_count << " gaps" << endl;

     // Write output.

     if (WRITE) 
     {    cout << Date( ) << ": building new unipaths" << endl;
          vecbasevector all(unibases);
          all.Append(merged_paths);
          for ( int id1 = 0; id1 < (int) unibases.size( ); id1++ ) 
          {    for (int j = 0; j < nexts[id1].isize( ); j++) 
               {    int id2 = nexts[id1][j];
	            basevector b = unibases[id1];
	            b.resize( b.size( ) + 1 );
	            b.Set( b.size( ) - 1, unibases[id2][K-1] );
	            all.push_back_reserve(b);    }    }
          vecKmerPath newpaths, newpathsrc, newunipaths;
          vec<tagged_rpint> newpathsdb, newunipathsdb;
          Mkdir777( run_dir + "/SimpleGapCloser.tmp" );
          ReadsToPathsCoreY( all, K, newpaths, newpathsrc, newpathsdb,
               run_dir + "/SimpleGapCloser.tmp", NUM_THREADS );
          Unipath( newpaths, newpathsrc, newpathsdb, newunipaths, newunipathsdb );
          KmerBaseBroker newkbb( K, newpaths, newpathsrc, newpathsdb, all );
          vecbasevector newunibases;
          for (size_t i = 0; i < newunipaths.size(); i++)
               newunibases.push_back_reserve( newkbb.Seq(newunipaths[i]) );
          cout << Date( ) << ": writing output files" << endl;
          String KS = ToString(K);
          String outhead = run_dir + "/" + HEAD_OUT;
          newpaths.WriteAll( outhead + ".paths.k" + KS );
          newpathsrc.WriteAll( outhead + ".paths_rc.k" + KS );
          newunipaths.WriteAll( outhead + ".unipaths.k" + KS );
          BinaryWrite3( outhead + ".pathsdb.k" + KS, newpathsdb );
          BinaryWrite3( outhead + ".unipathsdb.k" + KS, newunipathsdb );
          newunibases.WriteAll( outhead + ".unibases.k" + KS );    }
     cout << Date( ) << ": done" << endl;    }
