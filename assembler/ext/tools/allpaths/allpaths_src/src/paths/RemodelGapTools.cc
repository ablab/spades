///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>

#include "Basevector.h"
#include "CoreTools.h"
#include "PairsManager.h"
#include "ParallelVecUtilities.h"
#include "kmers/KmerRecord.h"
#include "math/Functions.h"
#include "paths/RemodelGapTools.h"
#include "random/NormalDistribution.h"

// GapComp
//
// d = possible gap
// a = unordered list of left co-starts
// b = unordered list of right stops
// x = ordered list of lengths
// X = histogram of lengths
// p = pairs (i,j) with i indexing an element of a, j indexing an element of b
// L1 = length of left contig
// L2 = length of right contig

vec<long double> GapComp( const vec<int> D, const vec<int>& a, const vec<int>& b, 
     const vec<int>& x, const vec<double>& X, const vec< pair<int,int> >& p, 
     const int L1, const int L2, const int VERBOSITY )
{    
     // Heuristics.

     const double xfloor = 0.0001;
     const long double no_data = -1000000000.0;

     // Computation.

     vec<long double> answer;
     int ma = Max(a), mb = Max(b), mx = X.isize( ) - 1;
     vec<int> A( ma + 1, 0 ), B( mb + 1, 0 );
     if ( VERBOSITY >= 2 ) PRINT(mx);
     for ( int i = 0; i < a.isize( ); i++ )
          A[ a[i] ]++;
     for ( int i = 0; i < b.isize( ); i++ )
          B[ b[i] ]++;
     int na = a.size( ), nb = b.size( ), nx = Sum(X);
     vec<int> mult( X.isize( ), 0 );
     for ( int i = 0; i < A.isize( ); i++ )
     {    if ( A[i] == 0 ) continue;
          for ( int j = 0; j < B.isize( ); j++ )
          {    if ( B[j] == 0 ) continue;
               if ( i + j < mult.isize( ) ) mult[i+j] += A[i] * B[j];    }    }
     for ( int q = 0; q < D.isize( ); q++ )
     {    int d = D[q];
          double N = 0.0;
          for ( int l = 0; l < X.isize( ); l++ )
          {    double xval;
               if ( l + d >= 0 && l + d < X.isize( ) ) xval = X[ l + d ];
               else xval = xfloor;
               N += mult[l] * xval / ( double(na)*double(nb)*double(nx) );    }
          if ( VERBOSITY >= 3 )
          {    cout << "\n";
               PRINT3( d, N, p.size( ) );    }
          double log_L = 0.0;
          vec<int> used1, used2;

          // Precheck for no data.

          int ndata = 0;
          for ( int m = 0; m < p.isize( ); m++ )
          {    int i = p[m].first, j = p[m].second;
               if ( a[i] + b[j] > mx ) continue;
               if ( a[i] + b[j] + d < mx ) ndata++;    }
          if ( ndata == 0 )
          {    answer.push_back(no_data);
               continue;    }

          // Proceed with computation.

          for ( int m = 0; m < p.isize( ); m++ )
          {    int i = p[m].first, j = p[m].second;
               if ( a[i] + b[j] > mx ) continue;
               used1.push_back(i), used2.push_back(j);
               double xval;
               int w = a[i] + b[j] + d;
               if ( w >= 0 && w < mx ) xval = X[w];
               else xval = xfloor;
               double add = log( A[ a[i] ] * B[ b[j] ] * xval
                    / ( N * double(na) * double(nb) * double(nx) ) );
               if ( VERBOSITY >= 4 ) PRINT3( a[i], b[j], add );
               log_L += log( A[ a[i] ] * B[ b[j] ] * xval
                    / ( N * double(na) * double(nb) * double(nx) ) );    }
          UniqueSort(used1), UniqueSort(used2);
          if ( VERBOSITY >= 3 ) PRINT(log_L);
          answer.push_back(log_L);    }
     return answer;    }

template<int K> void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          starts.push_back( starts.back( ) + Max( 0, u.isize( ) - K + 1 ) );    }
     kmers_plus.resize( starts.back( ) );
     for ( size_t i = 0; i < tigs.size( ); i++ )
     {    const basevector& u = tigs[i];
          #pragma omp parallel for
          for ( int jz = 0; jz <= u.isize( ) - K; jz += 1000 )
          {    kmer<K> x, xrc;
               for ( int j = jz; j <= Min( u.isize( ) - K, jz + 1000 ); j++ )
               {    int64_t r = starts[i] + j;
                    x.SetToSubOf( u, j ); 
                    xrc = x;
                    xrc.ReverseComplement( );
                    Bool fw = ( x < xrc );
                    kmers_plus[r].first = ( fw ? x : xrc );
                    kmers_plus[r].second = i; 
                    kmers_plus[r].third = ( fw ? j : -j-1 );    }    }    }
     ParallelSort(kmers_plus);    }

template<int K> void MakeKmerLookup0( const vecbasevector& unibases,
     vec< triple<kmer<K>,int,int> >& kmers_plus )
{    vec<int64_t> starts;
     starts.push_back(0);
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          starts.push_back( starts.back( ) + u.isize( ) - K + 1 );    }
     kmers_plus.resize( starts.back( ) );
     #pragma omp parallel for
     for ( size_t i = 0; i < unibases.size( ); i++ )
     {    const basevector& u = unibases[i];
          kmer<K> x;
          for ( int j = 0; j <= u.isize( ) - K; j++ )
          {    int64_t r = starts[i] + j;
               x.SetToSubOf( u, j ); 
               kmers_plus[r].first = x;
               kmers_plus[r].second = i; 
               kmers_plus[r].third = j;    }    }
     ParallelSort(kmers_plus);    }

template void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple< kmer<20>,int,int> >& kmers_plus );

template void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple< kmer<40>,int,int> >& kmers_plus );

template void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple< kmer<80>,int,int> >& kmers_plus );

template void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple< kmer<96>,int,int> >& kmers_plus );

template void MakeKmerLookup( const vecbasevector& tigs,
     vec< triple< kmer<100>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<20>,int,int> >& kmers_plus );

template void MakeKmerLookup0( const vecbasevector& tigs,
     vec< triple< kmer<96>,int,int> >& kmers_plus );

void FindTrueGaps( 
     // inputs:
     const vecbasevector& genome, const vecbasevector& tigs,
     const vec<superb>& scaffolds,
     // control:
     const vec<int>& tigs_to_process, const Bool TIME_STAMPS,
     // outputs:
     vec<int>& true_gap, vec<Bool>& true_gap_computed )
{
     int ngaps = 0;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
          ngaps += scaffolds[s].Ngaps( );
     true_gap.resize( ngaps, 0 );
     true_gap_computed.resize( ngaps, False );
     if (TIME_STAMPS) cout << Date( ) << ": computing gkmers_plus" << endl;
     const int L1 = 100;
     vec< triple<kmer<L1>,int,int> > gkmers_plus;
     MakeKmerLookup( genome, gkmers_plus );
     if (TIME_STAMPS) cout << Date( ) << ": forming gkmers" << endl;
     vec< kmer<L1> > gkmers( gkmers_plus.size( ) );
     for ( size_t i = 0; i < gkmers.size( ); i++ )
          gkmers[i] = gkmers_plus[i].first;
     if (TIME_STAMPS) cout << Date( ) << ": finding true gaps" << endl;
     int gap_id = -1;
     for ( int s = 0; s < scaffolds.isize( ); s++ )
     {    const superb& S = scaffolds[s];
          for ( int j = 0; j < S.Ngaps( ); j++ )
          {    gap_id++;
               int m1 = S.Tig(j), m2 = S.Tig(j+1);
               const basevector &M1 = tigs[m1], &M2 = tigs[m2];
               if ( tigs_to_process.nonempty( ) && !BinMember(tigs_to_process,m1) )
                    continue;
               kmer<L1> x1, x1rc, x2, x2rc;
               x1.SetToSubOf( M1, M1.isize( ) - L1 );
               x1rc = x1;
               x1rc.ReverseComplement( );
               Bool fw1 = ( x1 < x1rc );
               if ( !fw1 ) x1 = x1rc;
               int64_t low1 = LowerBound( gkmers, x1 );
               int64_t high1 = UpperBound( gkmers, x1 );
               x2.SetToSubOf( tigs[m2], 0 );
               x2rc = x2;
               x2rc.ReverseComplement( );
               Bool fw2 = ( x2 < x2rc );
               if ( !fw2 ) x2 = x2rc;
               int64_t low2 = LowerBound( gkmers, x2 );
               int64_t high2 = UpperBound( gkmers, x2 );
               for ( int pass = 1; pass <= 2; pass++ )
               {    if ( true_gap_computed[gap_id] ) continue;
                    const int L2 = ( pass == 1 ? 100 : 1000 );
                    if ( M1.isize( ) < L2 || M2.isize( ) < L2 ) continue;
                    vec<int64_t> hits1, hits2;
                    for ( int64_t l = low1; l < high1; l++ )
                    {    int pos1 = gkmers_plus[l].third; 
                         Bool xfw1 = ( pos1 >= 0 );
                         if ( !xfw1 ) pos1 = -pos1-1;
                         int g1 = gkmers_plus[l].second;
                         const basevector& G1 = genome[g1];
                         Bool mismatch = False;
                         if ( fw1 == xfw1 )
                         {    for ( int p = 0; p < L2 - L1; p++ )
                              {    int gpos1 = pos1 + p + L1 - L2;
                                   if ( gpos1 < 0 || M1[ M1.isize( ) - L2 + p ]
                                        != G1[gpos1] )
                                   {    mismatch = True;
                                        break;    }    }    }
                         else
                         {    for ( int p = 0; p < L2 - L1; p++ )
                              {    int gpos1 = pos1 + L2 - 1 - p;
                                   if ( gpos1 >= G1.isize( ) 
                                        || M1[ M1.isize( ) - L2 + p ]
                                        != 3 - G1[gpos1] )
                                   {    mismatch = True;
                                        break;    }    }    }
                         if ( !mismatch ) hits1.push_back(l);    }
                    if ( !hits1.solo( ) ) continue;
                    for ( int64_t l = low2; l < high2; l++ )
                    {    int pos2 = gkmers_plus[l].third; 
                         Bool xfw2 = ( pos2 >= 0 );
                         if ( !xfw2 ) pos2 = -pos2-1;
                         int g2 = gkmers_plus[l].second;
                         const basevector& G2 = genome[g2];
                         Bool mismatch = False;
                         if ( fw2 == xfw2 )
                         {    for ( int p = L1; p < L2; p++ )
                              {    int gpos2 = pos2 + p;
                                   if ( gpos2 >= G2.isize( ) 
                                        || M2[p] != G2[ gpos2] )
                                   {    mismatch = True;
                                        break;    }    }    }
                         else
                         {    for ( int p = L1; p < L2; p++ )
                              {    int gpos2 = pos2 + L1 - 1 - p;
                                   if ( gpos2 < 0 || M2[p] != 3 - G2[gpos2] )
                                   {    mismatch = True;
                                        break;    }    }    }
                         if ( !mismatch ) hits2.push_back(l);    }
                    if ( !hits2.solo( ) || gkmers_plus[ hits1[0] ].second 
                         != gkmers_plus[ hits2[0] ].second )
                    {    continue;    }
                    int pos1 = gkmers_plus[ hits1[0] ].third; 
                    int pos2 = gkmers_plus[ hits2[0] ].third;
                    Bool xfw1 = ( pos1 >= 0 ), xfw2 = ( pos2 >= 0 );
                    if ( !xfw1 ) pos1 = -pos1-1;
                    if ( !xfw2 ) pos2 = -pos2-1;
                    const int true_gap_upper_bound = 100000;
                    int g = ( fw1 == xfw1 ? pos2 - (pos1+L1) : pos1 - (pos2+L1) );
                    if ( ( int(fw1) + int(fw2) + int(xfw1) + int(xfw2) ) % 2 != 0
                         || Abs(g) > true_gap_upper_bound )
                    {    // PRINT8( m1, m2, int(fw1), int(fw2), 
                         // int(xfw1), int(xfw2), pos1, pos2 );
                         continue;    }
                    // cout << "gap from " << m1 << " to " 
                    //      << m2 << " is " << g << endl;
                    true_gap[gap_id] = g;
                    true_gap_computed[gap_id] = True;    }    }    }
     cout << Date( ) << ": found true size for " << Sum(true_gap_computed)
          << " gaps of " << true_gap_computed.size( ) << endl;    }

void PrintSideBySide( const String& z1, const String& z2, const int N )
{    istringstream zin1(z1), zin2(z2);
     String line1, line2;
     while(1)
     {    getline( zin1, line1 );
          getline( zin2, line2 );
          if ( !zin1.fail( ) && !zin2.fail( ) )
          {    cout << line1;
               for ( int i = 0; i < (int) N - (int) line1.size( ); i++ )
                    cout << " ";
               cout << line2 << "\n";    }
          else if ( zin1.fail( ) && !zin2.fail( ) )
          {    for ( int i = 0; i < N; i++ )
                    cout << " ";
               cout << line2 << "\n";    }
          else if ( !zin1.fail( ) && zin2.fail( ) ) cout << line1 << "\n";
          else break;    }    }

template<int K> void AlignReads( const Bool use_tail, const vecbasevector& tigs,
     const vecbasevector& bases, vec<Bool>& tried, vec<Bool>& placed_fw, 
     vec<Bool>& placed_rc, vec< pair<int,int> >& placement, const Bool TIME_STAMPS )
{
     if (TIME_STAMPS) cout << Date( ) << ": making kmer table" << endl;
     vec< triple<kmer<K>,int,int> > tig_kmers_plus;
     MakeKmerLookup( tigs, tig_kmers_plus );
     vec< kmer<K> > tig_kmers( tig_kmers_plus.size( ) );
     for ( size_t i = 0; i < tig_kmers.size( ); i++ )
          tig_kmers[i] = tig_kmers_plus[i].first;
     if (TIME_STAMPS) cout << Date( ) << ": forming read alignments" << endl;
     #pragma omp parallel for
     for ( size_t id = 0; id < bases.size( ); id++ )
     {    if ( tried[id] ) continue;
          if ( bases[id].isize( ) < K ) continue;
          tried[id] = True;
          kmer<K> x, xrc;
          int read_length = bases[id].size( );
          x.SetToSubOf( bases[id], ( use_tail ? read_length - K : 0 ) );
          int64_t low, high;
          Bool fw;
          GetBounds( tig_kmers, x, xrc, fw, low, high );
          if ( high - low != 1 ) continue;
          int tig = tig_kmers_plus[low].second, pos = tig_kmers_plus[low].third;
          Bool fw_tig = ( pos >= 0 );
          if ( !fw_tig ) pos = -pos-1;
          if ( fw == fw_tig ) placed_fw[id] = True;
          else placed_rc[id] = True;
          if (use_tail)
          {    if ( !fw_tig ) pos += K;
               if ( fw && fw_tig ) pos -= ( read_length - K );
               if ( fw && !fw_tig ) pos -= K;
               if ( !fw && !fw_tig ) pos -= read_length;    }
          else
          {    if ( fw != fw_tig ) pos -= ( read_length - K );    }
          placement[id].first = tig, placement[id].second = pos;    }    }

void AlignReads( const int K1, const int K2, const Bool use_tail, 
     const vecbasevector& tigs, const vecbasevector& bases, 
     const PairsManager& pairs, vec<Bool>& placed_fw, vec<Bool>& placed_rc, 
     vec< pair<int,int> >& placement, vec< vec<longlong> >& aligns_index, 
     const Bool TIME_STAMPS )
{
     vec<Bool> tried( bases.size( ), False );
     placed_fw.resize( bases.size( ), False );
     placed_rc.resize( bases.size( ), False );
     placement.resize( bases.size( ) );
     if ( K1 == 40 )
     {    AlignReads<40>( use_tail, tigs, bases, tried, placed_fw, placed_rc, 
               placement, TIME_STAMPS );    }
     else if ( K1 == 80 )
     {    AlignReads<80>( use_tail, tigs, bases, tried, placed_fw, placed_rc, 
               placement, TIME_STAMPS );    }
     else 
     {    cout << "AlignReads not implemented for K1 = " << K1 << endl;
          exit(1);    }
     Bool need2 = False;
     for ( size_t id = 0; id < bases.size( ); id++ )
          if ( bases[id].isize( ) >= K2 && bases[id].isize( ) < K1 ) need2 = True;
     if (need2) 
     {    if ( K2 == 20 )
          {    AlignReads<20>( use_tail, tigs, bases, tried, placed_fw, 
                    placed_rc, placement, TIME_STAMPS );    }
          else
          {    cout << "AlignReads not implemented for K2 = " << K2 << endl;
               exit(1);    }    }
     if (TIME_STAMPS) cout << Date( ) << ": indexing alignments" << endl;
     aligns_index.resize( tigs.size( ) );
     for ( size_t pid = 0; pid < pairs.nPairs( ); pid++ )
     {    int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
          if ( placed_fw[id1] || placed_rc[id1] )
               aligns_index[ placement[id1].first ].push_back(pid);
          if ( placement[id1].first == placement[id2].first ) continue; // Do not add the pair twice if aligned to same contig
          if ( placed_fw[id2] || placed_rc[id2] )
               aligns_index[ placement[id2].first ].push_back(pid);    }    }

template<int K> void PredictGap( 
     // inputs:
     const String& libtype,
     const vec< triple<kmer<K>,int,int> >& kmers_plus, const vec< kmer<K> >& kmers,
     const int max_overlap, const int nlibs, const int lib_offset,
     const vec<int>& libs_to_use, const vec<Bool>& lfail, const vecbasevector& tigs,
     const vec< vec<int> >& DISTS, const vec< vec<double> >& X, 
     const vecbasevector& bases, const PairsManager& pairs, 
     const vec<unsigned short>& read_lengths, 
     const vec< vec<longlong> >& aligns_index, const vec<Bool>& placed_fw, 
     const vec<Bool>& placed_rc, const vec< pair<int,int> >& placement, 
     const gap_id& gx, const int VERBOSITY, 
     // outputs:
     Bool& gap_predicted, int& gap, double& dev )
{
     const int M1 = ( gx.Type( ) == gap_id::BETWEEN ?  gx.M1() : gx.M() );
     const int M2 = ( gx.Type( ) == gap_id::BETWEEN ?  gx.M2() : gx.M() );
     vec<int> gaps;
     int step = 5;
     int mlow = -(max_overlap+step)/step * step;
     int top = 0;
     for ( int l = 0; l < nlibs; l++ ) 
          if (DISTS[l].size() > 0)
               top = Max( top, DISTS[l].back( ) );
     top = (top+step-1)/step * step;
     for ( int g = mlow; g <= top; g += step )
          gaps.push_back(g);
     vec< vec<long double> > ps;
     for ( int l = 0; l < nlibs; l++ )
     {    if ( libs_to_use.nonempty( ) && !Member( libs_to_use, l ) ) continue;
          if ( lfail[l] ) continue;
          int L1, L2;
          if ( gx.Type( ) == gap_id::BETWEEN )
          {    L1 = tigs[ gx.M1( ) ].size( ), L2 = tigs[ gx.M2( ) ].size( );    }
          else
          {    L1 = gx.Start( ), L2 = tigs[ gx.M( ) ].isize( ) - gx.Stop( );    }
          vec<int> a, b, x;
          vec< pair<int,int> > p;
          x = DISTS[l];
          vec<longlong> pids;
          if ( gx.Type( ) == gap_id::BETWEEN )
          {    pids = aligns_index[ gx.M1( ) ];
               pids.append( aligns_index[ gx.M2( ) ] );    }
          else pids = aligns_index[ gx.M( ) ];
          UniqueSort(pids);
          int unpairs = 0;
          for ( int pi = 0; pi < pids.isize( ); pi++ )
          {    longlong pid = pids[pi];
               int libid = lib_offset + pairs.libraryID(pid);
               if ( libid != l ) continue;
	       // Assign the pair id so that read(id1) is forward aligned on m1
	       // and read(id2) backward aligned on m2
               // int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
               int64_t id1 = -1, id2 = -1;
               for ( int pass = 1; pass <= 2; pass++ ) {
                    int64_t idx = ( pass == 1 ? pairs.ID1(pid): pairs.ID2(pid) );
                    if ( placement[idx].first == M1 && placed_fw[idx] )  id1 = idx;
                    if ( placement[idx].first == M2 && placed_rc[idx] )  id2 = idx;
               }
	       int d1 = ( id1 < 0 ? -1 : L1 - placement[id1].second );
	       int d2 = ( id2 < 0 ? -1 : placement[id2].second + read_lengths[id2] );
               if ( gx.Type( ) == gap_id::WITHIN ) d2 -= gx.Stop( );
               Bool use1, use2;
               if ( gx.Type( ) == gap_id::BETWEEN )
               {    use1 = ( id1 >=0 && placed_fw[id1] && placement[id1].first == gx.M1( )
                         && d1 <= x.back( ) );
                    use2 = ( id2 >=0 && placed_rc[id2] && placement[id2].first == gx.M2( )
                         && d2 <= x.back( ) );    }
               else
               {    use1 = ( id1 >=0 && placed_fw[id1] && placement[id1].first == gx.M( ) 
                         && d1 <= x.back( ) ) && d1 >= read_lengths[id1];
                    use2 = ( id2 >=0 && placed_rc[id2] && placement[id2].first == gx.M( )
                         && d2 <= x.back( ) ) && d2 >= read_lengths[id2];    }
               if ( use1 && use2 )
               {    p.push( a.size( ), b.size( ) );
                    a.push_back(d1);
                    b.push_back(d2);    }
               else
               {    const int max_errs = 5;
                    const int min_overlap = 40;
                    if (use1) 
                    {    a.push_back(d1);
                         if ( gx.Type( ) == gap_id::BETWEEN && libtype == "frag" )
                         {
                         if ( read_lengths[id2] >= K )
                         {    if ( !placed_rc[id2] 
                                   || placement[id2].first != gx.M1( ) ) 
                              {    kmer<K> x, xrc;
                                   x.SetToSubOf( bases[id2], 0 );
                                   xrc = x;
                                   xrc.ReverseComplement( );
                                   Bool fw2 = ( x < xrc );
                                   if ( !fw2 ) x = xrc;
                                   int64_t low = LowerBound( kmers, x );
                                   int64_t high = UpperBound( kmers, x );
                                   Bool found = False;
                                   int pos1 = placement[id1].second;
                                   for ( int64_t j = low; j < high; j++ )
                                   {    int pos2 = kmers_plus[j].third;
                                        Bool fw2_tig = ( pos2 >= 0 );
                                        if ( fw2 == fw2_tig ) continue;
                                        if ( !fw2_tig ) pos2 = -pos2 - 1 + K;
                                        if ( kmers_plus[j].second == gx.M1( ) ) 
                                        {    int d = pos2 - pos1;
                                             if ( d >= 0 && d <= DISTS[l].back( ) )
                                             {    found = True;
                                                  break;    }    }
                                        if ( kmers_plus[j].second == gx.M2( ) ) 
                                        {    int d = L1 - pos1 + pos2;
                                             // cout << "a: ";
                                             // PRINT5( id1, id2, pos1, pos2, d );
                                             if ( d >= 0 && d <= DISTS[l].back( ) )
                                             {    found = True;
                                                  break;    }    }    }
                                   if ( !found )
                                   {    const basevector& T1 = tigs[ gx.M1( ) ];
                                        for ( int stop1 = min_overlap; 
                                             stop1 < bases[id2].isize( ); stop1++ )
                                        {    int b2 = bases[id2].size( );
                                             int d = T1.isize( ) - pos1
                                                  + b2 - stop1;
                                             if ( d < 0 || d > DISTS[l].back( ) )
                                                  continue;
                                             int errs = 0;
                                             for ( int u = 0; u < stop1; u++ )
                                             {    if ( b2 - u - 1 < 0 ) continue;
                                                  if ( b2 - u - 1 >= b2 ) continue;
                                                  if ( u + T1.isize( ) - stop1 < 0 )
                                                       continue;
                                                  if ( u + T1.isize( ) - stop1
                                                       >= T1.isize( ) )
                                                  {    continue;    }
                                                  if ( 3 - bases[id2][ b2 - u  - 1 ]
                                                       !=
                                                       T1[ u + T1.isize( ) - stop1 ]
                                                       )
                                                  {    errs++;
                                                       if ( errs > max_errs )
                                                            break;    }    }
                                             if ( errs <= max_errs )
                                             {    found = True;
                                                  break;    }    }
                                        if ( !found )
                                        {    // PRINT2( id1, id2 );
                                             unpairs++;    }    }    }    }    }    }
                    if (use2) 
                    {    b.push_back(d2);    
                         if ( gx.Type( ) == gap_id::BETWEEN && libtype == "frag" )
                         {
                         if ( read_lengths[id1] >= K )
                         {    if ( !placed_fw[id1] 
                                   || placement[id1].first != gx.M2( ) )
                              {    kmer<K> x, xrc;
                                   x.SetToSubOf( bases[id1], 0 );
                                   xrc = x;
                                   xrc.ReverseComplement( );
                                   Bool fw1 = ( x < xrc );
                                   if ( !fw1 ) x = xrc;
                                   int64_t low = LowerBound( kmers, x );
                                   int64_t high = UpperBound( kmers, x );
                                   Bool found = False;
                                   int pos2 = placement[id2].second + K;
                                   for ( int64_t j = low; j < high; j++ )
                                   {    int pos1 = kmers_plus[j].third;
                                        Bool fw1_tig = ( pos1 >= 0 );
                                        if ( fw1 != fw1_tig ) continue;
                                        if ( !fw1_tig ) pos1 = -pos1 - 1 + K;
                                        if ( kmers_plus[j].second == gx.M2( ) ) 
                                        {    int d = pos2 - pos1;
                                             if ( d >= 0 && d <= DISTS[l].back( ) )
                                             {    found = True;
                                                  break;    }    }
                                        if ( kmers_plus[j].second == gx.M1( ) ) 
                                        {    int d = L1 - pos1 + pos2;
                                             // cout << "b: ";
                                             // PRINT5( id1, id2, pos1, pos2, d );
                                             if ( d >= 0 && d <= DISTS[l].back( ) )
                                             {    found = True;
                                                  break;    }    }    }
                                   if ( !found )
                                   {    const basevector& T2 = tigs[ gx.M2( ) ];
                                        for ( int start1 = 1; start1 < 
                                             bases[id1].isize( ) - min_overlap;
                                             start1++ )
                                        {    int d = pos2 + start1;
                                             if ( d < 0 || d > DISTS[l].back( ) )
                                                  continue;
                                             int errs = 0;
                                             for ( int u = start1; 
                                                  u < bases[id1].isize( ); u++ )
                                             {    if ( u-start1 >= T2.isize( ) )
                                                       break;
                                                  if ( bases[id1][u] !=
                                                       T2[u-start1] )
                                                  {    errs++;
                                                       if ( errs > max_errs )
                                                            break;    }    }
                                             if ( errs <= max_errs )
                                             {    found = True;
                                                  break;    }    }
                                        if ( !found )
                                        {    // PRINT2( id1, id2 );
                                             unpairs++;    
                                             }    }    }    }    }    }    }    }
          if ( VERBOSITY >= 2 ) PRINT4( a.size( ), b.size( ), p.size( ), unpairs );
          const int min_links = 1;
          if ( p.isize( ) < min_links ) continue;
          vec<long double> psl = GapComp(gaps, a, b, x, X[l], p, L1, L2, VERBOSITY);
          long double pmax = Max(psl), pssum = 0.0;
          for ( int i = 0; i < psl.isize( ); i++ )
          {    psl[i] = expl( psl[i] - pmax );
               pssum += psl[i];    }
          for ( int i = 0; i < psl.isize( ); i++ )
               psl[i] /= pssum;
          ps.push_back(psl);    }
     gap_predicted = ps.nonempty( );

     // Merge pdfs.

     vec<long double> PS( gaps.size( ), 1 );
     for ( int l = 0; l < ps.isize( ); l++ )
     {    for ( int j = 0; j < gaps.isize( ); j++ )
               PS[j] *= ps[l][j];    }
     long double sum = Sum(PS);
     for ( int j = 0; j < gaps.isize( ); j++ )
          PS[j] /= sum;

     // Compute gap and dev.  The "correct" computation of dev is commented out
     // because it doesn't seem to be better.

     gap = gaps[ FirstSumTo( PS, 0.5 ) ];
     double ndev = 3.0;
     double p = NormalDistribution( ndev, 0.0, 1.0 );
     int low = gaps[ FirstSumTo( PS, 1.0 - p ) ];
     int high = gaps[ FirstSumTo( PS, p ) ];
     dev = int( ceil( double(high-low) / ( 2.0 * ndev ) ) );
     /*
     double psum = 0.0;
     for ( int j = 0; j < gaps.isize( ); j++ )
          psum += PS[j] * double( gaps[j] - gap ) * double( gaps[j] - gap );
     dev = sqrt(psum);
     */

     // Multiple dev by a fudge factor which has the effect of raising low values:
     // dev   fudge
     //   5   2.22
     //  10   1.72
     // 100   1.25

     double fudge = 1.0 + 2.5 / pow( log( Max( dev, 2.0 )), 1.5 );
     dev = int( ceil( dev * fudge ) );

     // Print results.

     if ( VERBOSITY >= 2 ) cout << "\n";
     const double tail_toss = 0.0015;
     int low_print = 0, high_print = PS.size( ) - 1;
     double low_sum = 0.0, high_sum = 0.0;
     while( low_print < PS.isize( ) && low_sum < tail_toss )
     {    low_sum += PS[low_print];
          low_print++;    }
     low_print--;
     while( high_print >= 0 && high_sum < tail_toss )
     {    high_sum += PS[high_print];
          high_print--;    }
     high_print++;
     for ( int i = low_print; i <= high_print; i++ )
     {    int gap = gaps[i];
          double p = PS[i];
          if ( VERBOSITY >= 3 ) PRINT2( gap, p );    }    }

int FirstSumTo( const vec<long double>& x, long double bound )
{    int answer;
     long double total = 0.0;
     for ( int j = 0; j < x.isize( ); j++ )
     {    total += x[j];
          if ( total >= bound ) return j;    }
     return x.size( ) - 1;    }

void DefineDistributions(
     // inputs:
     const int nlibs, const int lib_offset, const vecbasevector& tigs, 
     const vec<unsigned short>& read_lengths, const PairsManager& pairs, 
     vec<Bool>& placed_fw, vec<Bool>& placed_rc, vec< pair<int,int> >& placement, 
     vec< vec<longlong> >& aligns_index, const Bool TIME_STAMPS,
     const Bool HISTOGRAMS,
     // outputs:
     vec<Bool>& lfail, vec<int>& max_dist, vec< vec<int> >& DISTS, 
     vec< vec<double> >& X )
{
     // Set a lower bound for contigs to use.  It's not clear that this makes sense.

     vec<int> len;
     int64_t total = 0, part = 0, min_to_use = 0;
     for ( size_t m = 0; m < tigs.size( ); m++ )
     {    len.push_back( tigs[m].size( ) );
          total += tigs[m].size( );    }
     ReverseSort(len);
     for ( size_t i = 0; i < len.size( ); i++ )
     {    part += len[i];
          min_to_use = len[i];
          if ( double(part)/double(total) >= 0.5 ) break;    }
     cout << Date( ) << ": "; PRINT(min_to_use);

     // Define visibility at given distances.

     int max_tig = 0;
     for ( size_t m = 0; m < tigs.size( ); m++ )
          max_tig = Max( max_tig, tigs[m].isize( ) );
     cout << Date( ) << ": maximum contig size = "
          << ToStringAddCommas(max_tig) << endl;
     vec<int64_t> VISIBLE( max_tig, 0 );
     for ( size_t m = 0; m < tigs.size( ); m++ )
     {    for ( size_t j = 0; j < tigs[m].size( ); j++ )
               VISIBLE[j] += tigs[m].isize( ) - j;    }

     // Compute distributions for libraries and upper cutoffs for them.  There are
     // two passes: first we find the cutoffs, then we compute the distributions,
     // taking into account the cutoffs.

     if (TIME_STAMPS) cout << Date( ) << ": computing library distributions" << endl;
     lfail.resize( nlibs, False );
     max_dist.resize(nlibs);
     DISTS.resize(nlibs);
     for ( int pass = 1; pass <= 2; pass++ )
     {    for ( int l = 0; l < nlibs; l++ )
               DISTS[l].clear( );
          for ( size_t m = 0; m < tigs.size( ); m++ )
          {    if ( tigs[m].isize( ) < min_to_use ) continue;
               int n = tigs[m].size( );
               vec<longlong> pids = aligns_index[m];
               UniqueSort(pids);
               for ( int pi = 0; pi < pids.isize( ); pi++ )
               {    longlong pid = pids[pi];
                    int64_t id1 = pairs.ID1(pid), id2 = pairs.ID2(pid);
                    if ( placement[id1].first != (int) m ) continue;
                    if ( placement[id2].first != (int) m ) continue;
                    if ( !placed_fw[id1] ) swap( id1, id2 );
                    if ( !placed_fw[id1] || !placed_rc[id2] ) continue;
                    int start1 = placement[id1].second;
                    int stop2 = placement[id2].second + read_lengths[id2];
                    if ( start1 < 0 || stop2 > n ) continue;
                    int dist = stop2 - start1;
                    if ( dist < 0 || dist >= max_tig ) continue;
                    int libid = lib_offset + pairs.libraryID(pid);
                    if ( pass == 2 && n - start1 < max_dist[libid] ) continue;
                    if ( pass == 2 && dist >= max_dist[libid] ) continue;
                    DISTS[libid].push_back(dist);    }    }
          for ( int l = 0; l < nlibs; l++ )
               Sort( DISTS[l] );

          if (HISTOGRAMS)
          {    if ( pass == 1 )
               {    for ( int l = 0; l < nlibs; l++ )
                    {    vec<int> bins(200, 0);
                         for ( int j = 0; j < DISTS[l].isize( ); j++ )
                         {    int x = DISTS[l][j]/100;
                              if ( x >= 200 ) break;
                              bins[x]++;    }
                         for ( int j = 0; j < bins.isize( ); j++ )
                              PRINT3( l, j, bins[j] );    }    }    }

          if ( pass == 2 ) break;

          // Define upper cutoff.  To do this, walk backwards from the end of the 
          // distribution, sampling 100 points at a time.  Taking into account
          // the visibility, estimate the noise rate and error bars for it, assuming
          // that the 100 points are all noise.  Now keep walking until we think
          // the signal:noise ratio is at least 10.  This defines the upper cutff.

          const double err_mult = 2.0;
          const double max_noise_frac = 0.1;
          const int bin_count = 100;
          for ( int l = 0; l < nlibs; l++ )
          {    lfail[l] = True;
               if ( DISTS[l].isize( ) >= bin_count )
               {    int N = DISTS[l].size( );
                    double min_freq_high = 1000000000;
                    for ( int rx = 0; N - (rx+1)*bin_count >= 0; rx++ )
                    {    double freq = 0, effective_size = 0;
                         int high = DISTS[l][ N - (rx)*bin_count - 1 ];
                         for ( int j = N - (rx+1)*bin_count; 
                              j < N - (rx)*bin_count; j++ )
                         {    int p = DISTS[l][j], q = high;
                              freq += Mean(VISIBLE) / double( VISIBLE[p] );
                              effective_size += 
                                   double(VISIBLE[q]) / double( VISIBLE[p] );    }
                         int low = DISTS[l][ N - (rx+1)*bin_count ];
                         freq /= double( high - low );
                         double freq_err = freq / sqrt(effective_size);
                         double freq_low = Max( 0.0, freq - err_mult * freq_err );
                         double freq_high = freq + err_mult * freq_err;
                         min_freq_high = Min( min_freq_high, freq_high );
                         double delta = freq_low / min_freq_high;
                         // PRINT5( low, high, freq_low, freq_high, delta );
                         if ( 1.0 / delta <= max_noise_frac )
                         {    max_dist[l] = high;
                              lfail[l] = False;
                              break;    }    }    }
               if ( lfail[l] && DISTS[l].nonempty( ) )
               {    int top = int( floor( double( DISTS[l].isize( ) * 0.99 ) ) );
                    max_dist[l] = DISTS[l][top];    
                    lfail[l] = False;    }    }    }
     vec<double> lib_mean(nlibs), lib_dev(nlibs);
     for ( int l = 0; l < nlibs; l++ )
     {    if ( !lfail[l] )
          {    lib_mean[l] = Mean( DISTS[l] );
               lib_dev[l] = StdDev( DISTS[l], lib_mean[l] );    }
          cout << Date( ) << ": library " << l << ", ";
          if ( lfail[l] ) cout << "can't find upper bound" << endl;
          else
          {    cout << lib_mean[l] << " +/- " << lib_dev[l] 
                    << ", upper bound = " << max_dist[l] << "\n";    }    }

     // Define distributions.

     const double xfloor = 0.0001;
     const double xfloor2 = 0.1;
     const int radius = 5;
     X.resize(nlibs);
     for ( int l = 0; l < nlibs; l++ )
     {    if ( lfail[l] ) continue;
          vec<int> x = DISTS[l];
          int mx = Max(x);
          X[l].resize( mx + 1, 0.0 );
          for ( int i = 0; i < x.isize( ); i++ )
          {    for ( int j = -radius; j <= radius; j++ )
               {    if ( x[i] + j >= 0 && x[i] + j < X[l].isize( ) )
                    {    X[l][ x[i] + j ]
                              += 1.0 / double( 2 * radius + 1 );    }    }    }
          for ( int i = 0; i < X[l].isize( ); i++ )
               X[l][i] = Max( xfloor2, X[l][i] );    }    }

void WriteReportAboutGaps( vec<int>& devs, vec<int>& adevs, vec<double>& offbys, 
     vec<double>& aoffbys, const vec< triple<int,int,int> >& results, 
     const vec<int>& true_gap, const vec<Bool>& true_gap_computed, 
     const int fails, const vec<int>& tigs_to_process )
{
     // Report count and devs.

     Sort(devs), Sort(adevs), Sort(offbys), Sort(aoffbys);
     int N = 25;
     int unassayed = 0;
     if ( tigs_to_process.empty( ) )
          unassayed = true_gap_computed.isize( ) - Sum(true_gap_computed);
     else
     {    for ( int j = 0; j < tigs_to_process.isize( ); j++ )
               if ( !true_gap_computed[ tigs_to_process[j] ] ) unassayed++;    }
     if ( devs.empty( ) )
     {    cout << "\ngaps:\n" << devs.size( ) << " + " << fails << " fails"
               << "\n+ " << unassayed << " unassayed\n";
          cout << "No gaps computed, report ends here." << endl;
          return;    }
     ostringstream outw1, outw2;
     outw1 << "\ngaps:\n" << devs.size( ) << " + " << fails << " fails"
           << "\n+ " << unassayed << " unassayed\n";
     int n_new = devs.size( ), n_old = adevs.size( );
     int new_25 = devs[ n_new / 4 ], old_25 = adevs[ n_old / 4 ];
     int new_50 = Median(devs), old_50 = Median(adevs);
     int new_75 = devs[ ( 3 * n_new ) / 4 ], old_75 = adevs[ ( 3 * n_old ) / 4 ];
     vec< vec<String> > rows;
     vec<String> row;
     row.push_back( "devs:", "new", "Q1 = " + ToString(new_25),
          "median = " + ToString(new_50), "Q3 = " + ToString(new_75) );
     rows.push_back(row);
     row.clear( );
     row.push_back( "", "old", "Q1 = " + ToString(old_25),
          "median = " + ToString(old_50), "Q3 = " + ToString(old_75) );
     rows.push_back(row);
     row.clear( );
     outw2 << "\n";
     PrintTabular( outw2, rows, 2, "lllll" );
     PrintSideBySide( outw1.str( ), outw2.str( ), N );
     double adj = double(Median(devs)) / double(Median(adevs));

     // Report off by in devs.  The last column old! is an adjusted
     // version, taking into account the difference between devs and adevs.

     cout << "\n";
     ostringstream outz1;
     outz1 << "percent off by:\n";
     rows.clear( );
     row.push_back( "dev", "new", "old", "old!" );
     rows.push_back(row);
     row.clear( );
     for ( double x = 2.0; x <= 5.0; x++ )
     {    int count = 0, acount = 0, aacount = 0;
          for ( int j = 0; j < offbys.isize( ); j++ )
          {    if ( offbys[j] > x ) count++;
               if ( aoffbys[j] > x ) acount++;
               if ( aoffbys[j] > x * adj ) aacount++;    }
          ostringstream out1, out2, out3;
          out1 << setiosflags(ios::fixed) << setprecision(1)
               << 100.0 * double(count) / double( offbys.size( ) );
          out2 << setiosflags(ios::fixed) << setprecision(1)
               << 100.0 * double(acount) / double( aoffbys.size( ) );
          out3 << setiosflags(ios::fixed) << setprecision(1)
               << 100.0 * double(aacount) / double( aoffbys.size( ) );
          row.push_back( ToString(x), out1.str( ), out2.str( ), out3.str( ) );  
          rows.push_back(row);
          row.clear( );    }
     PrintTabular( outz1, rows, 2, "rrrr" );
     String z1 = outz1.str( );

     // Report absolute off by.

     int err_lt1 = 0, err_lt10 = 0, err_lt100 = 0, err_lt1000 = 0, err_lt10000 = 0;
     int err_huge = 0;
     int xerr_lt1 = 0, xerr_lt10 = 0, xerr_lt100 = 0, xerr_lt1000 = 0;
     int xerr_lt10000 = 0, xerr_huge = 0;
     vec<int> errs;
     for ( int j = 0; j < results.isize( ); j++ )
     {    int err = Abs( results[j].first - results[j].second );
          errs.push_back( results[j].first - results[j].second );
          if ( err < 1 ) err_lt1++;
          else if ( err < 10 ) err_lt10++;
          else if ( err < 100 ) err_lt100++;
          else if ( err < 1000 ) err_lt1000++;
          else if ( err < 10000 ) err_lt10000++;
          else err_huge++;
          int xerr = Abs( results[j].third - results[j].second );
          if ( xerr < 1 ) xerr_lt1++;
          else if ( xerr < 10 ) xerr_lt10++;
          else if ( xerr < 100 ) xerr_lt100++;
          else if ( xerr < 1000 ) xerr_lt1000++;
          else if ( xerr < 10000 ) xerr_lt10000++;
          else xerr_huge++;    }
     rows.clear( );
     row.push_back( "error", "new", "old" );
     rows.push_back(row), row.clear( );
     row.push_back( "< 1", ToString(err_lt1), ToString(xerr_lt1) );
     rows.push_back(row), row.clear( );
     row.push_back( "< 10", ToString(err_lt10), ToString(xerr_lt10) );
     rows.push_back(row), row.clear( );
     row.push_back( "< 100", ToString(err_lt100), ToString(xerr_lt100) );
     rows.push_back(row), row.clear( );
     row.push_back( "< 1000", ToString(err_lt1000), ToString(xerr_lt1000) );
     rows.push_back(row), row.clear( );
     row.push_back( "< 10000", ToString(err_lt10000), ToString(xerr_lt10000) );
     rows.push_back(row), row.clear( );
     row.push_back( ">= 10000", ToString(err_huge), ToString(xerr_huge) );
     rows.push_back(row), row.clear( );
     ostringstream outz2;
     outz2 << "absolute error:\n";
     PrintTabular( outz2, rows, 2, "lrr" );
     String z2 = outz2.str( );

     // Merge tables.

     PrintSideBySide( z1, z2, N );
     Sort(errs);
     cout << "\nmedian absolute error = " << Median(errs) << endl;    }

void GetLibCount( const String& run_dir, const String& head, 
     int& nlibs, const int VERBOSITY, const Bool TIME_STAMPS )
{
     if (TIME_STAMPS) cout << Date( ) << ": loading " << head << endl;
     String pairs_file = run_dir + "/" + head;
     if ( head == "frag_reads_filt" ) pairs_file += "_cpd.pairs";
     else pairs_file += ".pairs";
     longlong nreads;
     vec<String> lib_names;
     vec<int> lib_sep, lib_sd;
     ReadPairsManagerLibInfo( pairs_file, nreads, lib_names, lib_sep, lib_sd );
     UniqueSort(lib_names);
     int nlibs_this = lib_names.size( );
     if ( VERBOSITY >= 1 )
     {    for ( int l = 0; l < nlibs_this; l++ )
          {    cout << Date( ) << ": library " << nlibs + l << " is of type "
                    << head.Before( "_reads" ) + "_reads" << "\n";    }    }
     nlibs += nlibs_this;    }

void GetLibCounts( const String& run_dir, const vec<String>& lib_types_to_use,
		   const int VERBOSITY, const Bool TIME_STAMPS,
		   const String& frags, const String& jumps, const String& long_jumps, 
		   int& nlibs_frag, int& nlibs_jump, int& nlibs_long )
{    int nlibs = 0;
     if ( Member( lib_types_to_use, String("frag") ) )
          GetLibCount( run_dir, frags, nlibs, VERBOSITY, TIME_STAMPS );
     nlibs_frag = nlibs;
     if ( Member( lib_types_to_use, String("jump") ) )
          GetLibCount( run_dir, jumps, nlibs, VERBOSITY, TIME_STAMPS );
     nlibs_jump = nlibs - nlibs_frag;
     if ( Member( lib_types_to_use, String("long") ) )
     {    if ( IsRegularFile( run_dir + "/" + long_jumps + ".pairs" ) )
	 {    GetLibCount( run_dir, long_jumps, 
                    nlibs, VERBOSITY, TIME_STAMPS );    }    }
     nlibs_long = nlibs - nlibs_frag - nlibs_jump;    }

template<int K> void ComputeOverlaps( const basevector& M1, 
     const basevector& M2, const int m1, const vec<int>& max_dist, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus, const vec< kmer<K> >& kmers,
     const int VERBOSITY, vec<int>& accepted_overlaps, int& max_overlap )
{
     // Test for plausibility of negative gaps.  Accept overlap o if at least
     // (o-K+1)/2 kmer hits are found within 10% of o.  A danger of this is that if
     // a contig has a hanging end which "overlaps" the other contig, we don't
     // report a negative overlap.  Another problem is that this won't see contigs
     // that are out of order.
     //
     // More serious problem: gaps near the floor have compressed standard
     // deviation.

     const double req_over_frac = 0.5;
     const double over_delta = 0.1;
     const int overlap_fudge = 40;
     vec<int> overlaps;
     for ( int p2 = 0; p2 <= Min( Max(max_dist), M2.isize( ) ) - K; p2++ )
     {    kmer<K> x, xrc;
          x.SetToSubOf( M2, p2 );
          xrc = x;
          xrc.ReverseComplement( );
          Bool fw2 = ( x < xrc );
          if ( !fw2 ) x = xrc;
          int64_t low = LowerBound( kmers, x ), high = UpperBound( kmers, x );
          for ( int64_t l = low; l < high; l++ )
          {    int m = kmers_plus[l].second;
               if ( m != m1 ) continue;
               int p1 = kmers_plus[l].third;
               Bool fw1 = ( p1 >= 0 );
               if ( fw1 != fw2 ) continue;
               if ( !fw1 ) p1 = -p1-1;
               overlaps.push_back( p2 + ( M1.isize( ) - p1 ) );    }    }
     Sort(overlaps);
     for ( int j = 0; j < overlaps.isize( ); j++ )
     {    int o = overlaps[j];
          int support = 1;
          for ( int l = j - 1; l >= 0; l-- )
          {    if ( o - overlaps[l] > double(o) * over_delta ) break;
               support++;    }
          for ( int l = j + 1; l < overlaps.isize( ); l++ )
          {    if ( overlaps[l] - o > double(o) * over_delta ) break;
               support++;    }
          if ( support >= (o-K+1) * req_over_frac ) accepted_overlaps.push_back(o);
          int k = overlaps.NextDiff(j);
          j = k - 1;    }
     max_overlap = overlap_fudge;
     if ( accepted_overlaps.nonempty( ) ) max_overlap += accepted_overlaps.back( );
     if ( VERBOSITY >= 2 ) 
     {    cout << "accepted overlaps:\n";
          for ( int j = 0; j < accepted_overlaps.isize( ); j++ )
               cout << accepted_overlaps[j] << "\n";    }    }

template void ComputeOverlaps( const basevector& M1,
     const basevector& M2, const int m1, const vec<int>& max_dist,
     const vec< triple<kmer<40>,int,int> >& kmers_plus, const vec< kmer<40> >& kmers,
     const int VERBOSITY, vec<int>& accepted_overlaps, int& max_overlap );

template void ComputeOverlaps( const basevector& M1,
     const basevector& M2, const int m1, const vec<int>& max_dist,
     const vec< triple<kmer<80>,int,int> >& kmers_plus, const vec< kmer<80> >& kmers,
     const int VERBOSITY, vec<int>& accepted_overlaps, int& max_overlap );

template void PredictGap( 
     const String& libtype,
     const vec< triple<kmer<40>,int,int> >& kmers_plus, const vec< kmer<40> >& kmers,
     const int max_overlap, const int nlibs, const int lib_offset,
     const vec<int>& libs_to_use, const vec<Bool>& lfail, const vecbasevector& tigs,
     const vec< vec<int> >& DISTS, const vec< vec<double> >& X, 
     const vecbasevector& bases,
     const PairsManager& pairs, const vec<unsigned short>& read_lengths, 
     const vec< vec<longlong> >& aligns_index, const vec<Bool>& placed_fw, 
     const vec<Bool>& placed_rc, const vec< pair<int,int> >& placement, 
     const gap_id& gx, const int VERBOSITY, 
     Bool& gap_predicted, int& gap, double& dev );

template void PredictGap( 
     const String& libtype,
     const vec< triple<kmer<80>,int,int> >& kmers_plus, const vec< kmer<80> >& kmers,
     const int max_overlap, const int nlibs, const int lib_offset,
     const vec<int>& libs_to_use, const vec<Bool>& lfail, const vecbasevector& tigs,
     const vec< vec<int> >& DISTS, const vec< vec<double> >& X, 
     const vecbasevector& bases,
     const PairsManager& pairs, const vec<unsigned short>& read_lengths, 
     const vec< vec<longlong> >& aligns_index, const vec<Bool>& placed_fw, 
     const vec<Bool>& placed_rc, const vec< pair<int,int> >& placement, 
     const gap_id& gx, const int VERBOSITY, 
     Bool& gap_predicted, int& gap, double& dev );
