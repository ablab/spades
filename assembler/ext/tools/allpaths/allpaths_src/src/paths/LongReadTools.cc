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
#include "FastIfstream.h"
#include "PackAlign.h"
#include "PrintAlignment.h"
#include "Set.h"
#include "graph/DigraphTemplate.h"
#include "kmers/KmerRecord.h"
#include "pairwise_aligners/SmithWatBandedA.h"
#include "pairwise_aligners/SmithWatFree.h"
#include "paths/AssemblyEdit.h"
#include "paths/LongReadTools.h"

int DistX(const GapPatcher0& p1, 
          const GapPatcher0& p2,
          const basevector& b1, 
          const basevector& b2,
          alignment * p_al,
          double * timer)
{
  // Extend p1 and p2 to the left and right so that they start and
  // end at the same places.
  basevector bp1 = p1.r, bp2 = p2.r;
  if      (p1.upos1 < p2.upos1) bp2 = Cat(basevector(b1, p1.upos1, p2.upos1 - p1.upos1), p2.r);
  else if (p1.upos1 > p2.upos1) bp1 = Cat(basevector(b1, p2.upos1, p1.upos1 - p2.upos1), p1.r);
  if      (p1.upos2 > p2.upos2) bp2 = Cat(bp2, basevector(b2, p2.upos2, p1.upos2 - p2.upos2));
  else if (p1.upos2 < p2.upos2) bp1 = Cat(bp1, basevector(b2, p1.upos2, p2.upos2 - p1.upos2));
  
  // Align them.
  int best_loc;
  int mismatch_penalty = 1;
  int gap_penalty = 1;
  const int score = (bp1.size() <= bp2.size()) ?
    SmithWatFree(bp1, bp2, best_loc, *p_al, true, true, mismatch_penalty, gap_penalty) :
    SmithWatFree(bp2, bp1, best_loc, *p_al, true, true, mismatch_penalty, gap_penalty);

  return score;
}

void DeleteX(GapPatcher0 * p_p, 
             int * p_sum, 
             vec<alignment> * p_als,
             const vec<GapPatcher0>& P, 
             const int L,
             const basevector& b1, 
             const basevector& b2,
             double * timer)
{
  basevector X(p_p->r, L, p_p->r.isize() - 2*L);
  for (int j = 0; j < X.isize(); j++) {
    cout << "delete j = " << j << " 0";

    basevector x = X;
    x.resize(j);
    x = Cat(x, basevector(X, j+1, X.isize() - (j+1)));
    GapPatcher0 px(*p_p);
    px.r = Cat(basevector(p_p->r, 0, L), 
               x, 
               basevector(p_p->r, p_p->r.isize() - L, L));

    int sumx = 0;
    #pragma omp parallel for
    for (int i = 0; i < P.isize(); i++) {
      alignment a;
      int d = DistX(px, P[i], b1, b2, &a, timer);
      #pragma omp critical
      { sumx += d; }
    }
    if (sumx < *p_sum) {
      *p_sum = sumx;
      X = x;    
      j--;
      cout << " 1";
    }
    else cout << " 0";
    cout << "   szx = " << X.size() 
         << "  *p_sum = " << *p_sum
         << "  sumx = " << sumx 
         << "  dsum = " << (sumx - *p_sum) 
         << endl; 
  }
  p_p->r = Cat(basevector(p_p->r, 0, L), 
               X, 
               basevector(p_p->r, p_p->r.isize() - L, L));
}

void InsertX(GapPatcher0 * p_p, 
             int * p_sum, 
             vec<alignment> * p_als, 
             const vec<GapPatcher0>& P, 
             const int L,
             const basevector& b1, 
             const basevector& b2,
             double * timer)
{
  basevector X(p_p->r, L, p_p->r.isize() - 2*L);
  for (int j = 0; j <= X.isize(); j++) {
    for (int k = 0; k < 4; k++) {
      cout << "insert j = " << j << " " << k;

      basevector x = X;
      x.resize(j);
      x.push_back(k);
      x = Cat(x, basevector(X, j, X.isize() - j));
      GapPatcher0 px(*p_p);
      px.r = Cat(basevector(p_p->r, 0, L), 
                 x, 
                 basevector(p_p->r, p_p->r.isize() - L, L));

      int sumx = 0;
      #pragma omp parallel for
      for (int i = 0; i < P.isize(); i++) {
        alignment a;
        int d = DistX(px, P[i], b1, b2, &a, timer);
        #pragma omp critical
        { sumx += d; }
      }
      if (sumx < *p_sum) {
        *p_sum = sumx;
        X = x;
        cout << " 1";
      }
      else cout << " 0";

      cout << "   szx = " << X.size() 
           << "  *p_sum = " << *p_sum
           << "  sumx = " << sumx  
           << "  dsum = " << (sumx - *p_sum) 
           << endl; 
    }
  }
  p_p->r = Cat(basevector(p_p->r, 0, L), 
               X, 
               basevector(p_p->r, p_p->r.isize() - L, L));
}

void SubstituteX(GapPatcher0 * p_p, 
                 int * p_sum, 
                 vec<alignment> * p_als, 
                 const vec<GapPatcher0>& P, 
                 const int L,
                 const basevector& b1, 
                 const basevector& b2,
                 double * timer)
{
  basevector X(p_p->r, L, p_p->r.isize() - 2*L);
  for (int j = 0; j < X.isize(); j++) {
    for (int k = 1; k < 4; k++) {
      cout << "substitute j = " << j << " " << k;

      basevector x = X;
      x.Set(j, (x[j] + k) % 4);
      GapPatcher0 px(*p_p);
      px.r = Cat(basevector(p_p->r, 0, L), 
                 x,
                 basevector(p_p->r, p_p->r.isize() - L, L));

      int sumx = 0;
      #pragma omp parallel for
      for (int i = 0; i < P.isize(); i++) {
        alignment a;
        int d = DistX(px, P[i], b1, b2, &a, timer);
        #pragma omp critical
        { sumx += d; }
      }
      if (sumx < *p_sum) {
        *p_sum = sumx;
        X = x;
        cout << " 1";
      }
      else cout << " 0";

      cout << "   szx = " << X.size() 
           << "  *p_sum = " << *p_sum
           << "  sumx = " << sumx  
           << "  dsum = " << (sumx - *p_sum) 
           << endl; 
    }
  }
  p_p->r = Cat(basevector(p_p->r, 0, L), 
               X, 
               basevector(p_p->r, p_p->r.isize() - L, L));
}

GapPatcher0 patcher_optimal_original2(const vec<GapPatcher0> & p0s,
                                      const size_t i_best, 
                                      const int L,
                                      const BaseVec & bv1,
                                      const BaseVec & bv2,
                                      double * timer) 
{
  const size_t np = p0s.size();

  // ---- Find minimum upos1 and maximum upos2
  int upos1_min = p0s[0].upos1;
  int upos2_max = p0s[0].upos2;
  for (size_t ip = 0; ip != np; ip++) {
    if (p0s[ip].upos1 < upos1_min) upos1_min = p0s[ip].upos1;
    if (p0s[ip].upos2 > upos2_max) upos2_max = p0s[ip].upos2;
  }
 

  // ---- Extend all BaseVecs to upos1_min and upos2_max
  vec<GapPatcher0> ps;
  for (size_t ip = 0; ip != np; ip++) {
    BaseVec bv = Cat(BaseVec(bv1, upos1_min, p0s[ip].upos1 - upos1_min), 
                     p0s[ip].r,
                     BaseVec(bv2, p0s[ip].upos2, upos2_max - p0s[ip].upos2));
    ps.push_back(GapPatcher0(bv, upos1_min, upos2_max));
  }



  GapPatcher0 X = ps[i_best];
  vec<alignment> as(np);
  // Compute the distance of p1 to the pack
  int sum = 0;
  for (size_t ip = 0; ip < np; ip++)
    sum += DistX(X, ps[ip], bv1, bv2, &(as[ip]), timer);

  // Optimize the distance
  const int max_reps = 10;
  for ( int r = 0; r < max_reps; r++ ) {
    cout << "pass= " << r << endl;
    GapPatcher0 X_save = X;
    DeleteX     (&X, &sum, &as, ps, L, bv1, bv2, timer);
    InsertX     (&X, &sum, &as, ps, L, bv1, bv2, timer);
    SubstituteX (&X, &sum, &as, ps, L, bv1, bv2, timer);
    if (X == X_save) break;
    X_save = X;
  }
  return X;
}



GapPatcher0 patcher_optimal_original(const vec<GapPatcher0> & p0s,
                                     const size_t i_best, 
                                     const int L,
                                     const BaseVec & bv1,
                                     const BaseVec & bv2,
                                     double * timer) 
{
  const size_t np = p0s.size();

  GapPatcher0 X = p0s[i_best];
  vec<alignment> as(np);
  // Compute the distance of p1 to the pack
  int sum = 0;
  for (size_t ip = 0; ip < np; ip++)
    sum += DistX(X, p0s[ip], bv1, bv2, &(as[ip]), timer);

  // Optimize the distance
  const int max_reps = 10;
  for ( int r = 0; r < max_reps; r++ ) {
    cout << "pass= " << r << endl;
    GapPatcher0 X_save = X;
    SubstituteX (&X, &sum, &as, p0s, L, bv1, bv2, timer);
    DeleteX     (&X, &sum, &as, p0s, L, bv1, bv2, timer);
    InsertX     (&X, &sum, &as, p0s, L, bv1, bv2, timer);
    if (X == X_save) break;
    X_save = X;
  }
  return X;
}







int KmerId(const basevector& b, const int L, const int p)
{
  int n = 0;
  for (int l = 0; l < L; l++) {
    n <<= 2;  // same as *= 4
    n += b[p + l];
  }
  return n;
}

Bool PerfectMatch( const basevector& b1, const basevector& b2,
     const int p1, const int p2, const int len )
{    for ( int i = 0; i < len; i++ )
     {    if ( p1+i >= b1.isize( ) || p2+i >= b2.isize( ) ) return False;
          if ( b1[p1+i] != b2[p2+i] ) return False;    }
     return True;    }

void GetLocalAligns( const basevector& r, const vecbasevector& U, 
     const vec< vec< pair<int,int> > >& Ulocs, const int L, const int flank,
     const int max_errs, vec< triple<int,int,int> >& aligns, ostream& out, 
     const Bool SHOW_ALIGNS, const Bool require_extend )
{
     basevector b, c;
     alignment a;
     for ( int p = 0; p <= r.isize( ) - L; p++ )
     {    int n = KmerId( r, L, p );
          for ( int v = 0; v < Ulocs[n].isize( ); v++ )
          {    int i = Ulocs[n][v].first;
               int j = Ulocs[n][v].second;
               if (require_extend)
               {    if ( j - p >= 0 && j - p + r.isize( ) <= U[i].isize( ) ) 
                    {    continue;    }    }

               // Now in principle we compare the region from p - flank to 
               // p + L + flank on the read to the region from j - flank to 
               // j + L + flank on the unibase.  If either of these goes 
               // "off the end", we shift flanking bases from one end to the 
               // other.  If this isn't possible, we give up.

               int left_flank = Min( p, j, flank );
               int right_flank 
                    = Min( r.isize( ) - L - p, U[i].isize( ) - L - j, flank );
               if ( left_flank < flank ) right_flank += flank - left_flank;
               if ( right_flank < flank ) left_flank += flank - right_flank;
               if ( p < left_flank || j < left_flank ) continue;
               if ( r.isize( ) - L - p < right_flank ) continue;
               if ( U[i].isize( ) - L - j < right_flank ) continue;
               b.SetToSubOf( r, p - left_flank, L + left_flank + right_flank );
               int best_loc;
               int errs = -1;

               // The region on the unibase to which we align is chosen to
               // be slightly larger than the region on the read.

               const int extra = 10;
               int left_extra = Min( extra, j - left_flank );
               int right_extra 
                    = Min( extra, U[i].isize( ) - L - j - right_flank );

               c.SetToSubOf( U[i], j - left_flank - left_extra, 
                    L + left_flank + right_flank + left_extra + right_extra );
               errs = SmithWatFree( b, c, best_loc, a, false, false, 1, 1 );
               errs = ActualErrors( b, c, a, 1, 1 );
               if ( errs < 0 || errs > max_errs ) continue;

               if (SHOW_ALIGNS)
               {    out << "i = " << i << ", j = " << j << ", p = " << p
                         << ", j-p = " << j-p 
                         << ", U[i].size = " << U[i].size( ) << ", b = ";
                    for ( int l = 0; l < b.isize( ); l++ )
                    {    if ( l == flank || l == flank+L ) out << " ";
                         out << as_base( b[l] );    }
                    out << ", errs = " << errs << "\n";    }

               aligns.push( i, p, j-p );    }    }
     // Sort(aligns);    
          }

void GetGlobalAligns( const basevector& r, const vecbasevector& U, 
     const vec< triple<int,int,int> >& aligns, const int bandwidth_div,
     vec<align>& aligns_a, const double sub_frac, const double ins_frac,
     const double del_frac )
{
     aligns_a.resize( aligns.size( ) );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    int u = aligns[i].first;
          int rpos = aligns[i].second, offset = aligns[i].third;
          int predicted_overlap = IntervalOverlap( 
               0, U[u].isize( ), offset, offset + r.isize( ) );
          int bandwidth = predicted_overlap / bandwidth_div;
          int errors;
          int sub = int(round(1.0/sub_frac));
          int del = int(round(1.0/ins_frac));
          int ins = int(round(1.0/del_frac));
          // int sub = 1;  /* 19? */
          // int ins = 1;  /* 45? */
          // int del = 1;  /* 8?  */
          sub = 1;
          ins = 1;
          del = 1;
          SmithWatBandedA2<unsigned int>( U[u], r, offset, bandwidth, 
               aligns_a[i], errors, (ostream*) 0, sub, ins, del );    }    }

void MakeAlignmentGraph( const vec< pair< int, vec< pair<int,int> > > >& alignsx,
     const vecbasevector& U, const int verbosity, ostream& rout, digraphE<int>& G )
{
     int N = alignsx.size( );
     vec< vec<int> > from(N), to(N), from_edge_obj(N), to_edge_obj(N);
     vec<int> overlaps;
     vec< triple<int,int,int> > ALIGNS;

     vec< triple<int,int,int> > flat;
     for ( int i = 0; i < alignsx.isize( ); i++ )
     {    for ( int j = 0; j < alignsx[i].second.isize( ); j++ )
          {    flat.push( alignsx[i].second[j].second,
                    alignsx[i].second[j].first, i );    }    }
     Sort(flat);
     vec< vec< pair<int,int> > > offsets(N);
     for ( int i = 0; i < flat.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < flat.isize( ); j++ )
               if ( flat[j].first != flat[i].first ) break;
          for ( int k1 = i; k1 < j; k1++ )
          for ( int k2 = i; k2 < j; k2++ )
          {    if ( k2 == k1 ) continue;
               int j1 = flat[k1].third, j2 = flat[k2].third;
               int v1 = alignsx[j1].first, v2 = alignsx[j2].first;
               int o = flat[k1].second - flat[k2].second;
               if ( U[v1].isize( ) >= U[v2].isize( ) + o ) continue;
               if ( v1 == v2 && o == 0 ) continue;
               Bool match = True;
               for ( int x1 = Max( 0, o ); x1 < U[v1].isize( ); x1++ )
               {    int x2 = x1 - o;
                    if ( x2 >= U[v2].isize( ) ) break;
                    if ( U[v1][x1] != U[v2][x2] )
                    {    match = False;
                         break;    }    }
               if (match) offsets[j1].push( j2, o );    }
          i = j - 1;    }
     for ( int j1 = 0; j1 < N; j1++ )
     {    UniqueSort( offsets[j1] );
          for ( int l = 0; l < offsets[j1].isize( ); l++ )
          {    int o = offsets[j1][l].second, j2 = offsets[j1][l].first;
               from[j1].push_back(j2), to[j2].push_back(j1);
               from_edge_obj[j1].push_back( overlaps.size( ) );
               to_edge_obj[j2].push_back( overlaps.size( ) );
               overlaps.push_back(o);
               int v1 = alignsx[j1].first, v2 = alignsx[j2].first;
               if ( verbosity >= 1 ) 
                    PRINT5_TO( rout, j1, j2, v1, v2, o );    }    }
     G.Initialize( from, to, overlaps, to_edge_obj, from_edge_obj );    }

void AlignToRef( const vecbasevector& B, const String& run_dir,
     const String& data_dir, ostream& rout, const String& module )
{    temp_file fastb( run_dir + "/tmp/" + module + ".XXXXXX" );
     B.WriteAll(fastb);
     fast_pipe_ifstream in( "QueryLookupTable K=12 MM=12 MC=0.15 "
          "SEQS=" + fastb + " SEQS_IS_FASTB=True L=" + data_dir + "/genome.lookup "
          "VISUAL=True SMITH_WAT=True NH=True QUIET=True" );
     String line;
     while(1)
     {    getline( in, line );
          if ( in.fail( ) ) break;
          rout << line << "\n";    }
     flush(rout);    }

// ConvertPathIntoBases
// p = a sequence of edges in graph G
// alignsx = as documented in LongReadTools.h for MakeAlignmentGraph

void ConvertPathIntoBases( const int L, const digraphE<int>& G, const vec<int>& p, 
     vec< pair< int, vec< pair<int,int> > > >& alignsx, const vecbasevector& U, 
     triple< basevector, int, int >& B )
{    vec<int> to_left, to_right;
     G.ToLeft(to_left), G.ToRight(to_right);
     const pair< int, vec< pair<int,int> > >& a1 = alignsx[ to_left[ p.front( ) ] ];
     const pair< int, vec< pair<int,int> > >& a2 = alignsx[ to_right[ p.back( ) ] ];
     basevector b = U[a1.first];
     for ( int j = 0; j < p.isize( ); j++ )
     {    int offset = G.EdgeObject( p[j] );
          int u = alignsx[ to_right[ p[j] ] ].first;
          int uprev = alignsx[ to_left[ p[j] ] ].first;
          b.resize( b.isize( ) - ( U[uprev].isize( ) - offset ) );
          b = Cat( b, U[u] );
          if ( j == p.isize( ) - 1 )
          {    int ustart = a1.second.front( ).first;
               int ustop = a2.second.back( ).first + L;

               // what's happening here?
               // p[j] = last edge in p
               // goes between two alignments
               // u = unipath in second alignment
               // ustop = uposn for second alignment + L
               // ustart = upos1 for second alignment
               // (again see alignsx documentation)

               // Sanity check.  Not that the purpose of 
               // ConvertPathIntoBasesValid is to precheck this.

               ForceAssertLt( ustart + U[u].isize( ) - ustop, b.isize( ) );

               B.first.SetToSubOf( 
                    b, ustart, b.isize( ) - ustart - ( U[u].isize( ) - ustop ) );
               B.second = a1.second.front( ).second;
               B.third = a2.second.back( ).second + L;    }    }    }

Bool ConvertPathIntoBasesValid( const int L, const digraphE<int>& G, 
     const vec<int>& p, vec< pair< int, vec< pair<int,int> > > >& alignsx, 
     const vecbasevector& U )
{    vec<int> to_left, to_right;
     G.ToLeft(to_left), G.ToRight(to_right);
     const pair< int, vec< pair<int,int> > >& a1 = alignsx[ to_left[ p.front( ) ] ];
     const pair< int, vec< pair<int,int> > >& a2 = alignsx[ to_right[ p.back( ) ] ];
     basevector b = U[a1.first];
     for ( int j = 0; j < p.isize( ); j++ )
     {    int offset = G.EdgeObject( p[j] );
          int u = alignsx[ to_right[ p[j] ] ].first;
          int uprev = alignsx[ to_left[ p[j] ] ].first;
          b.resize( b.isize( ) - ( U[uprev].isize( ) - offset ) );
          b = Cat( b, U[u] );
          if ( j == p.isize( ) - 1 )
          {    int ustart = a1.second.front( ).first;
               int ustop = a2.second.back( ).first + L;
               return ( ustart + U[u].isize( ) - ustop < b.isize( ) );    }    }
     return True;    }

void ConvertPathsIntoBases( const int L, const basevector& target,
     const digraphE<int>& G, const vec< vec<int> >& epaths,
     const vec<int>& vpaths, vec< pair< int, vec< pair<int,int> > > >& alignsx, 
     const vecbasevector& U, vec<basevector>& bpaths, const int verbosity, 
     ostream& rout, const String& run_dir, const String& data_dir,
     const String& module, const Bool simple )
{
     vec<int> to_left, to_right;
     G.ToLeft(to_left), G.ToRight(to_right);
     vec< triple< basevector, int, int > > Bs;
     for ( int i = 0; i < epaths.isize( ); i++ )
     {    const vec<int>& p = epaths[i];
          triple< basevector, int, int > B;
          ConvertPathIntoBases( L, G, p, alignsx, U, B );
          Bs.push_back(B);    }
     for ( int i = 0; i < vpaths.isize( ); i++ )
     {    triple< basevector, int, int > B;
          const pair< int, vec< pair<int,int> > >& a = alignsx[ vpaths[i] ];
          int u = a.first;
          int ustart = a.second.front( ).first, ustop = a.second.back( ).first + L;
          B.first.SetToSubOf( U[u], ustart, ustop - ustart );
          B.second = a.second.front( ).second;
          B.third = a.second.back( ).second + L;
          Bs.push_back(B);    }
     if (simple) 
     {    for ( int j = 0; j < Bs.isize( ); j++ )
               bpaths.push_back( Bs[j].first );
          return;    }
     UniqueSort(Bs);
     vec< vec<int> > from( Bs.size( ) ), to( Bs.size( ) );
     for ( int i1 = 0; i1 < Bs.isize( ); i1++ )
     {    for ( int i2 = 0; i2 < Bs.isize( ); i2++ )
          {    if ( Bs[i1].third <= Bs[i2].second )
               {    from[i1].push_back(i2);
                    to[i2].push_back(i1);    }    }    }
     digraph H( from, to );
     vec< vec<int> > hpaths;
     H.AllPaths( -1, -1, hpaths );
     vec<int> len( hpaths.size( ), 0 );
     for ( int j = 0; j < hpaths.isize( ); j++ )
     {    for ( int l = 0; l < hpaths[j].isize( ); l++ )
               len[j] += Bs[ hpaths[j][l] ].third - Bs[ hpaths[j][l] ].second;    }
     ReverseSortSync( hpaths, len );
     for ( int j = 0; j < hpaths.isize( ); j++ )
     {    if ( j > 0 && len[j] < len[j-1] ) break;
          basevector b(target);
          const vec<int>& h = hpaths[j];
          for ( int l = h.isize( ) - 1; l >= 0; l-- )
          {    basevector b1( b, 0, Bs[ h[l] ].second );
               basevector b2( Bs[ h[l] ].first );
               basevector b3( b, Bs[ h[l] ].third, -1 );
               b = Cat( b1, b2, b3 );    }
          bpaths.push_back(b);    }
     UniqueSort(bpaths);
     if ( verbosity >= 1 ) PRINT_TO( rout, bpaths.size( ) );
     if ( verbosity >= 3 && bpaths.nonempty( ) )
     {    rout << "\nalignments of bpaths:\n";
          vecbasevector B( bpaths.size( ) );
          for ( int j = 0; j < bpaths.isize( ); j++ )
               B[j] = bpaths[j];
          #pragma omp critical
          {    AlignToRef( B, run_dir, data_dir, rout, module );    }    }    }

int SmithWatFreeSym( const basevector& b1, const basevector& b2, align& a,
     const Bool penalize_left_gap, const Bool penalize_right_gap,
     unsigned int mismatch_penalty, unsigned int gap_penalty )
{    
     alignment al;
     int best_loc, errs;
     if ( b1.size( ) <= b2.size( ) )
     {    errs = SmithWatFree( b1, b2, best_loc, al, 
               penalize_left_gap, penalize_right_gap,
               mismatch_penalty, gap_penalty );
          a.UnpackFrom(al);    }
     else
     {    errs = SmithWatFree( b2, b1, best_loc, al, 
               penalize_left_gap, penalize_right_gap,
               mismatch_penalty, gap_penalty );
          a.UnpackFrom(al);
          a.Flip( );    }
     return errs;    }

vec<placementx> FindGenomicPlacements( basevector b, const int L,
     const vecbasevector& genome, const vec< vec< pair<int,int> > >& Glocs )
{    vec<placementx> places;
     for ( int pass = 1; pass <= 2; pass++ )
     {    if ( pass == 2 ) b.ReverseComplement( );
          int n = KmerId( b, L, 0 );
          for ( int z = 0; z < Glocs[n].isize( ); z++ )
          {    const basevector& g = genome[ Glocs[n][z].first ];
               const int gpos = Glocs[n][z].second;
               if ( !PerfectMatch( b, g, 0, gpos, b.size( ) ) ) continue;
               places.push( Glocs[n][z].first, gpos, pass == 1 );    }    }
     return places;    }





void patchers_gap_collect(const BaseVecVec & bvs,
                          const vec<GapPatcher> & patchers, 
                          const unsigned np_min, 
                          const int sz_padding_min,
                          vec<vec<GapPatcher> > * p_patchers_gap)
{
  const unsigned np = patchers.size();
  vec<GapPatcher> patchers_tmp;
  for (unsigned ip = 0; ip != np; ip++) {
    
    // ---- only include patchers that leave enough padding space on 
    //      the reads they are patching 

    if (patchers[ip].tpos1 >= sz_padding_min &&
        bvs[patchers[ip].t2].isize() - patchers[ip].tpos2 >= sz_padding_min)
      patchers_tmp.push_back(patchers[ip]);

    const bool last_in_group = (ip == np - 1 ||
                                (patchers[ip].t1 != patchers[ip+1].t1 ||
                                 patchers[ip].t2 != patchers[ip+1].t2));
    if (last_in_group) {
      if (patchers_tmp.size() >= np_min)
        p_patchers_gap->push_back(patchers_tmp);
      patchers_tmp.clear();
    }
  }
}




// type == 0   choose best patcher based on shortest patcher size 
// type == 1   choose best patcher based on median patcher size 
// type == 2   choose best patcher based on median gap size 

void patchers_first_guesses(const BaseVecVec & T,
                            const vec<vec<GapPatcher> > & patchers_gap,
                            const unsigned type,
                            const unsigned sz_patcher_min, // only for type 0
                            vec<unsigned> * p_ips_best)
{
  const unsigned n_gaps = patchers_gap.size();
  for (unsigned ig = 0; ig != n_gaps; ig++) {
    const unsigned np = patchers_gap[ig].size();
    
    if (type == 2) {   // choose best patcher based on median gap size 
      const int nb1 = T[patchers_gap[ig][0].t1].size();
      vec<unsigned> ips;
      vec<int> gaps;
      for (unsigned ip = 0; ip != np; ip++) {
        const GapPatcher & p = patchers_gap[ig][ip];
        gaps.push_back(int(p.r.size()) - p.tpos2 + p.tpos1 - nb1);
        ips.push_back(ip);
      }
      SortSync(gaps, ips);
      p_ips_best->push_back(ips[np/2]);  // pick the index of the median 
    }
    
    else if (type == 1) {  // choose best patcher based on median patcher size 
      const int nb1 = T[patchers_gap[ig][0].t1].size();
      vec<unsigned> ips;
      vec<int> szs;
      for (unsigned ip = 0; ip != np; ip++) {
        const GapPatcher & p = patchers_gap[ig][ip];
        szs.push_back(p.r.size());
        ips.push_back(ip);
      }
      SortSync(szs, ips);
      p_ips_best->push_back(ips[np/2]);  // pick the index of the median 
    }
    
    else if (type == 0) {  // choose best patcher based on shortest patcher size 
      unsigned sz_min = 0;
      unsigned ip_shortest = np;  // np is an invalid ip
      for (unsigned ip = 0; ip != np; ip++) {
        const GapPatcher & p = patchers_gap[ig][ip];
        if ((p.r.size() >= sz_patcher_min) && 
            (ip_shortest == np ||    // ip_shortest still invalid?          
             p.r.size() < sz_min)) {
          sz_min = p.r.size();
          ip_shortest = ip;
        }
      }   
      p_ips_best->push_back(ip_shortest);
    }
  }
  
}                                               




// ---- sort gaps according to their predicted performance
//      time scales as  np . len^2 
// 
void indexes_gaps_for_parallel(const vec<vec<GapPatcher> >& patchers_gap, 
                               const vec<unsigned> & ips_best,
                               vec<unsigned> * p_is_gaps,
                               const bool verbose)
{
  const unsigned n_gaps = patchers_gap.size();
  vec<double> time;
  for (unsigned ig = 0; ig != n_gaps; ig++) {
    const double np = patchers_gap[ig].size();
    const size_t ip = ips_best[ig];
    if (ip < np) {  // make sure ip is valid
      const double len = patchers_gap[ig][ip].r.size();
      time.push_back(np * len * len);
    }
    else { // invalid ip will be done last (and discarded)
      time.push_back(0.0);
    }
    p_is_gaps->push_back(ig);
  }
  ReverseSortSync(time, *p_is_gaps); // longer patches are done last
  if (verbose) {
    for (unsigned ig = 0; ig != n_gaps; ig++) {
      cout << "ig_sorted= " << setw(3) << ig
           << " ig= " << (*p_is_gaps)[ig] 
           << " time= " << setw(10) << time[ig] 
           << endl;
    }
  }
}
     


void patcher_short_print(const int i_gap,
                         const GapPatcher0 & p_opt, 
                         const BaseVec & bv1, 
                         const BaseVec & bv2)
{
  uint nb1 = bv1.size();
  uint nb2 = bv2.size();
  uint nbp = p_opt.r.size();
  uint upos1 = p_opt.upos1;
  uint upos2 = p_opt.upos2;
  int gap = (int(nbp) + int(upos1) - int(nb1) - int(upos2)); 
  cout  << "i_gap= " << i_gap
        << "  upos1= " << setw(6) << upos1 
        << "  bv1.size()= " << setw(8) << nb1 
        << "  upos2= " << setw(6) << upos2 
        << "  bv2.size()= " << setw(8) << nb2 
        << "  r.size()= " << setw(8) << nbp 
        << "  gap= " << setw(6) << gap
        << endl;
             
  uint nb1_left  = Min(10u, upos1);
  uint nb1_right = Min(20u, nb1 - upos1);

  uint nb2_left  = Min(20u, upos2);
  uint nb2_right = Min(10u, nb2 - upos2);
             
  cout << "i_gap= " << i_gap << " ... ";
  for (uint i = 0; i < nb1_left; i++) cout << as_base(bv1[upos1 - nb1_left + i]);
  cout << " ";
  for (uint i = 0; i < nb1_right; i++) cout << as_base(bv1[upos1 + i]);
  cout << " ... gap = " << setw(4) << gap << " ... ";
  for (uint i = 0; i < nb2_left; i++) cout << as_base(bv2[upos2 - nb2_left + i]);
  cout << " ";
  for (uint i = 0; i < nb2_right; i++) cout << as_base(bv2[upos2 + i]);
  cout << endl;

  
  uint nbp_right = Min(nb1_right, nbp);
  uint nbp_left  = Min(nb2_left,  nbp);

  cout << "i_gap= " << i_gap << "     ";
  for (uint i = 0; i < nb1_left; i++) cout << " ";
  cout << " ";
  for (uint i = 0; i < nbp_right; i++) cout << as_base(p_opt.r[i]);
  for (uint i = nbp_right; i < nb1_right; i++) cout << "_";
  cout << " .................. ";
  for (uint i = nbp_left; i < nb2_left; i++) cout << "_";
  for (uint i = 0; i < nbp_left; i++) cout << as_base(p_opt.r[nbp - nbp_left + i]);
  cout << " ";
  for (uint i = 0; i < nb1_right; i++) cout << " ";
  cout << endl;
}

template<int K> 
void CleanPatch(

     // inputs:

     const int L, const vecbasevector& U, 
     const vec< vec< pair<int,int> > >& Ulocs, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus,
     const fastavector& LEFT, 
     const fastavector& RIGHT,

     // inputs and outputs:

     assembly_edit& e,

     // logging:

     ostringstream& outx, 
     const Bool DIRECT, 
     const int verbosity, 
     const String& data_dir, 
     const String& run_dir )
{
     int start1 = e.Start1( ), stop2 = e.Stop2( );
     ForceAssertEq( e.Nreps( ), 1 );
     basevector patch = e.Rep(0);
     fastavector left, right;
     left.SetToSubOf( LEFT, 0, start1 );
     right.SetToSubOf( RIGHT, stop2, (int) RIGHT.size( ) - stop2 );

     if ( (int) left.size( ) < K || (int) right.size( ) < K )
     {    if ( verbosity >= 1 ) outx << "flanks too short" << endl;
          return;    }

     // Find flanking kmers.  Define the target.

     vec<flank_datum> F;
     GetFlankData( left, right, patch, kmers_plus, F, verbosity, outx, 
          run_dir, data_dir );

     if ( F.empty( ) )
     {    if ( verbosity >= 1 ) 
               outx << "failed to find flanking kmer, giving up" << endl;
          return;    }    

     // Go through the flank data.  Note that the way we use ambiguous bases to 
     // define multiple flanks is probably not the best way to do things.  Rather,
     // when we make the original long read patches, we should use the long reads
     // to determine these ambiguous bases.

     vec<basevector> best_bpaths;
     Bool exploded = False;
     if ( verbosity >= 1 )
          outx << Date( ) << ": found " << F.size( ) << " flanks" << endl;
     for ( int fi = 0; fi < F.isize( ); fi++ )
     {    const flank_datum& f = F[fi];

          // Define the alignments corresponding to the end unipaths.

          vec< pair< int, vec< pair<int,int> > > > alignsx;
          vec<double> mismatch_ratesx;
          vec< pair<int,int> > M;
          basevector u1_trunc( U[f.u1], f.pos1, U[f.u1].isize( ) - f.pos1 );
          align a;
          vec<ho_interval> pf1, pf2;
          SmithWatFreeSym( u1_trunc, f.target, a, True, False );
          a.PerfectIntervals1( u1_trunc, f.target, pf1 );
          a.PerfectIntervals2( u1_trunc, f.target, pf2 );
          for ( int j = 0; j < pf1.isize( ); j++ )
          {    int start1 = pf1[j].Start( ), stop1 = pf1[j].Stop( );
               int start2 = pf2[j].Start( );
               for ( int l = start1; l <= stop1 - L; l++ )
                    M.push( f.pos1 + l, start2 + l - start1 );    }
          alignsx.push( f.u1, M );
          mismatch_ratesx.push_back(0);
          M.clear( );
          basevector u2_trunc( U[f.u2], 0, f.pos2 + K );
          SmithWatFreeSym( u2_trunc, f.target, a, False, True );
          a.PerfectIntervals1( u2_trunc, f.target, pf1 );
          a.PerfectIntervals2( u2_trunc, f.target, pf2 );
          for ( int j = 0; j < pf1.isize( ); j++ )
          {    int start1 = pf1[j].Start( ), stop1 = pf1[j].Stop( );
               int start2 = pf2[j].Start( );
               for ( int l = start1; l <= stop1 - L; l++ )
                    M.push( l, start2 + l - start1 );    }
          alignsx.push( f.u2, M );
          mismatch_ratesx.push_back(0);

          // Search for L-mer matches to the patch.

          int START, STOP;
          if ( verbosity >= 1 )
          {    outx << Date( ) << ": flank " << fi+1 
                    << ", aligning to target" << endl;    }
          AlignToTarget( f.target, K, L, U, Ulocs, alignsx, mismatch_ratesx, 
               START, STOP, verbosity, outx );

          // Form graph and find paths, then convert paths into bases, then
          // pick the best path.

          digraphE<int> G;
          if ( verbosity >= 1 )
          {    outx << Date( ) << ": flank " << fi+1 
                    << ", making alignment graph" << endl;    }
          MakeAlignmentGraph( alignsx, U, verbosity, outx, G );
          vec< vec<int> > epaths;
          vec<int> vpaths;
          const int max_paths = 100000;
          const int max_iterations = 100000;
          if ( verbosity >= 1 ) outx << Date( ) << ": creating edgepaths" << endl;
          Bool OK = G.EdgePaths(START, STOP, epaths, -1, max_paths, max_iterations);
          if ( !OK )
          {    if ( verbosity >= 1 ) 
                    outx << "edge path computation exploded" << endl;
               break;    }
          if ( verbosity >= 1 ) 
               outx << Date( ) << ": found " << epaths.size( ) << " paths" << endl;
          Bool exploded = False;
          if ( epaths.empty( ) )
          {    if ( verbosity >= 1 ) 
               {    outx << "no complete paths, looking for "
                         << "incomplete ones" << endl;    }
               for ( int v = 0; v < G.N( ); v++ )
                    if ( G.Source(v) && G.Sink(v) ) vpaths.push_back(v);
               for ( int v = 0; v < G.N( ); v++ )
               {    if (exploded) break;
                    if ( !G.Source(v) ) continue;
                    for ( int w = 0; w < G.N( ); w++ )
                    {    if ( !G.Sink(w) ) continue;
                         if ( v == START && w == STOP ) continue;
                         vec< vec<int> > epathsi;
                         Bool OK = G.EdgePaths( 
                              v, w, epathsi, -1, max_paths, max_iterations );
                         if ( !OK )
                         {    if ( verbosity >= 1 ) 
                              {    outx << "edge path computation exploded"
                                        << endl;    }
                              exploded = True;
                              break;    }

                         // Remove nonsensical paths.

                         vec<Bool> to_remove( epathsi.size( ), False );
                         for ( int j = 0; j < epathsi.isize( ); j++ )
                         {    if ( !ConvertPathIntoBasesValid( L, G, epathsi[j],
                                   alignsx, U ) )
                              {    to_remove[j] = True;    }    }
                         EraseIf( epathsi, to_remove );

                         // Save paths.

                         epaths.append(epathsi);    }    }    }
          if ( !exploded )
          {    if ( verbosity >= 1 ) 
                    PRINT2_TO( outx, epaths.size( ), vpaths.size( ) );
               vec<basevector> bpaths;
               ConvertPathsIntoBases( L, f.target, G, epaths, vpaths, alignsx, 
                    U, bpaths, verbosity, outx, run_dir, data_dir,
                    "CleanLongReadPatches" );
               if ( bpaths.nonempty( ) ) 
               {    basevector bpath;
                    PickBestPath( bpaths, f.target, bpath, verbosity, 
                         outx, run_dir, data_dir );
                    best_bpaths.push_back(bpath);    }    }    }

     // Update the edit.  For now we require that the best bpaths all have
     // the same size and have no 'new' ambiguous bases.

     if ( !exploded && best_bpaths.nonempty( ) ) 
     {    start1 -= F[0].lshift + K, stop2 += F[0].rshift + K;
          Bool same_size = True;
          int n = best_bpaths[0].size( );
          for ( int j = 1; j < best_bpaths.isize( ); j++ )
               if ( best_bpaths[j].isize( ) != n ) same_size = False;
          if (same_size)
          {    fastavector X( best_bpaths[0] );
               for ( int i = 1; i < best_bpaths.isize( ); i++ )
                    X.combine( fastavector( best_bpaths[i] ) );
               while ( X.size( ) > 0 && start1 < (int) LEFT.size( ) - 1 )
               {    if ( LEFT[start1] != X[0] ) break;
                    X.SetToSubOf( X, 1, X.size( ) - 1 );
                    start1++;    }
               while( X.size( ) > 0 && stop2 > 0 )
               {    if ( RIGHT[stop2-1] != X[ X.size( )-1 ] ) break;
                    X.resize( X.size( ) - 1 );
                    stop2--;    }
               if ( X.AmbCount( ) == 0 )
               {    if ( verbosity >= 1 ) outx << "saving patch" << endl;
                    e.SetStart1(start1), e.SetStop2(stop2);
                    e.Rep(0) = X.ToBasevector( );    }    }    }    }

template void CleanPatch<96>(
     const int L, const vecbasevector& U, const vec< vec< pair<int,int> > >& Ulocs, 
     const vec< triple<kmer<96>,int,int> >& kmers_plus, const fastavector& LEFT, 
     const fastavector& RIGHT, assembly_edit& e, ostringstream& outx, 
     const Bool DIRECT, const int verbosity, const String& data_dir,
     const String& run_dir );

template void CleanPatch<640>(
     const int L, const vecbasevector& U, const vec< vec< pair<int,int> > >& Ulocs, 
     const vec< triple<kmer<640>,int,int> >& kmers_plus, const fastavector& LEFT, 
     const fastavector& RIGHT, assembly_edit& e, ostringstream& outx, 
     const Bool DIRECT, const int verbosity, const String& data_dir,
     const String& run_dir );

template< int K > void GetFlankData( const fastavector& left, 
     const fastavector& right, const basevector& patch, 
     const vec< triple<kmer<K>,int,int> >& kmers_plus, 
     vec<flank_datum>& F, const int verbosity, ostream& rout, 
     const String& run_dir, const String& data_dir )
{
     const int max_amb = 2;
     flank_datum f;
     kmer<K> x;
     fastavector Lx, Rx;
     vec< triple<int,int,basevector> > left_data, right_data;
     for ( f.lshift = 0; f.lshift < (int) left.size( ) - K + 1; f.lshift++ )
     {    Lx.SetToSubOf( left, (int) left.size( ) - K - f.lshift, K );
          if ( Lx.AmbCount( ) > max_amb ) return;
          vecbasevector LEFTS = Lx.AllBasevectors( );
          for ( size_t j = 0; j < LEFTS.size( ); j++ )
          {    x.Set( LEFTS[j] );
               int64_t p = BinPosition1( kmers_plus, x );
               if ( p >= 0 ) 
                    left_data.push( kmers_plus[p].second, kmers_plus[p].third,
                         LEFTS[j] );    }
          if ( left_data.nonempty( ) ) break;    }
     for ( f.rshift = 0; f.rshift < (int) right.size( ) - K + 1; f.rshift++ )
     {    Rx.SetToSubOf( right, f.rshift, K );
          if ( Rx.AmbCount( ) > max_amb ) return;
          vecbasevector RIGHTS = Rx.AllBasevectors( );
          for ( size_t j = 0; j < RIGHTS.size( ); j++ )
          {    x.Set( RIGHTS[j] );
               int64_t p = BinPosition1( kmers_plus, x );
               if ( p >= 0 ) 
                    right_data.push( kmers_plus[p].second, kmers_plus[p].third,
                         RIGHTS[j] );    }
          if ( right_data.nonempty( ) ) break;    }
     for ( int j1 = 0; j1 < left_data.isize( ); j1++ )
     for ( int j2 = 0; j2 < right_data.isize( ); j2++ )
     {    f.u1 = left_data[j1].first, f.pos1 = left_data[j1].second;
          f.u2 = right_data[j2].first, f.pos2 = right_data[j2].second;
          if ( verbosity >= 1 ) 
          {    rout << Date( ) << ": ";
               PRINT2_TO( rout, f.u1, f.u2 );    }
          basevector L0, R0;
          fastavector L0x, R0x;
          L0x.SetToSubOf( left, (int) left.size( ) - K - f.lshift, K + f.lshift );
          R0x.SetToSubOf( right, 0, K + f.rshift );
          L0 = L0x.ToBasevector( ), R0 = R0x.ToBasevector( );
          f.target = Cat( L0, patch, R0 );
          for ( int i = 0; i < K; i++ )
          {    f.target.Set( i, left_data[j1].third[i] );
               f.target.Set( i + f.target.isize( ) - K, 
                    right_data[j2].third[i] );    }
          F.push_back(f);
          if ( verbosity >= 2 )
          {    rout << "\nalignment of target:\n";
               vecbasevector T(1);
               T[0] = f.target;
               #pragma omp critical
               {    AlignToRef( T, run_dir, data_dir, 
                         rout, "CleanLongReadPatches" );    }    }    }    }

void AlignToTarget( 
     /* inputs */
     const basevector& target, const int K, const int L, const vecbasevector& U,
     const vec< vec< pair<int,int> > >& Ulocs,
     /* outputs */
     vec< pair< int, vec< pair<int,int> > > >& alignsx, 
     vec<double>& mismatch_ratesx, int& START, int& STOP,
     /* logging */
     const int verbosity, ostream& rout 
          )
{
     vec< triple<int,int,int> > aligns;
     vec<double> mismatch_rates;
     alignment al;
     int best_loc;
     set< triple<int,int,int> > Mused;
     const int max_freq = 10000;
     for ( int p = 0; p <= target.isize( ) - L; p++ )
     {    int n = KmerId( target, L, p );
          if ( Ulocs[n].isize( ) > max_freq ) continue;
          for ( int v = 0; v < Ulocs[n].isize( ); v++ )
          {    int u = Ulocs[n][v].first, up = Ulocs[n][v].second;
               if ( Member( Mused, make_triple(u,up,p) ) ) continue;
               basevector uleft( U[u], 0, up );
               basevector uright( U[u], up + L, U[u].isize( ) - up - L );
               basevector tleft( target, 0, p );
               basevector tright( target, p + L, target.isize( ) - p - L );
               if ( uleft.size( ) > tleft.size( ) ) continue;
               if ( uright.size( ) > tright.size( ) ) continue;
               int left_bound, right_bound, mis1 = 0, mis2 = 0;
               vec<ho_interval> perf1, perf2, PERF1, PERF2;
               if ( uleft.size( ) > 0 )
               {    SmithWatFree( uleft, tleft, best_loc, al, False, True );
                    align a(al);
                    left_bound = a.pos2( );
                    vec<int> mgg = a.MutationsGap1Gap2( uleft, tleft );
                    mis1 = mgg[0];
                    a.PerfectIntervals1( uleft, tleft, PERF1 );
                    a.PerfectIntervals2( uleft, tleft, PERF2 );    }
               else left_bound = p;
               if ( PERF1.nonempty( ) && PERF1.back( ).Stop( ) == up
                    && PERF2.back( ).Stop( ) == p )
               {    PERF1.back( ).AddToStop(L), PERF2.back( ).AddToStop(L);    }
               else
               {    PERF1.push( up, up + L ), PERF2.push( p, p + L );    }
               if ( uright.size( ) > 0 )
               {    SmithWatFree( uright, tright, best_loc, al, True, False );
                    align a(al);
                    right_bound = a.Pos2( ) + tleft.isize( ) + L;    
                    vec<int> mgg = a.MutationsGap1Gap2( uright, tright );
                    mis2 = mgg[0];
                    a.PerfectIntervals1( uright, tright, perf1 );
                    a.PerfectIntervals2( uright, tright, perf2 );
                    for ( int j = 0; j < perf1.isize( ); j++ )
                    {    perf1[j].Shift( uleft.isize( ) + L );
                         perf2[j].Shift( tleft.isize( ) + L );    }
                    if ( perf1.nonempty( ) 
                         && PERF1.back( ).Stop( ) == perf1[0].Start( )
                         && PERF2.back( ).Stop( ) == perf2[0].Start( ) )
                    {    PERF1.back( ).SetStop( perf1[0].Stop( ) );
                         PERF2.back( ).SetStop( perf2[0].Stop( ) );
                         perf1.erase( perf1.begin( ) );
                         perf2.erase( perf2.begin( ) );    }
                    PERF1.append(perf1), PERF2.append(perf2);    }
               else right_bound = tleft.isize( ) + L;
               const double max_mismatch_rate = 0.1;
               double mismatch_rate = double( mis1 + mis2 ) / double( U[u].size( ) );
               if ( mismatch_rate > max_mismatch_rate ) continue;
               vec< pair<int,int> > M;
               for ( int j = 0; j < PERF1.isize( ); j++ )
               {    int start1 = PERF1[j].Start( ), stop1 = PERF1[j].Stop( );
                    int start2 = PERF2[j].Start( );
                    for ( int l = start1; l <= stop1 - L; l++ )
                         M.push( l, start2 + l - start1 );    }
               alignsx.push( u, M );
               mismatch_ratesx.push_back(mismatch_rate);
               const int max_offset_diff = 5;
               for ( int j = 0; j < M.isize( ); j++ )
               {    int offset_diff = Abs( up-p - (M[j].first-M[j].second) );
                    if ( offset_diff <= max_offset_diff )
                    {    Mused.insert( make_triple(
                              u, M[j].first, M[j].second ) );    }    }
               aligns.push( left_bound, right_bound, u );    
               mismatch_rates.push_back(mismatch_rate);    }    }
     UniqueSortSync( aligns, mismatch_rates );
     if ( verbosity >= 1 )
     {    for ( int j = 0; j < aligns.isize( ); j++ )
          {    rout << "unipath " << aligns[j].third << " aligned "
                    << aligns[j].first << "-" << aligns[j].second << " of " 
                    << target.size( ) << ", error rate = " << setprecision(3)
                    << 100.0 * mismatch_rates[j] << "% " << endl;    }    }

     // Merge alignments.

     equiv_rel E( alignsx.size( ) );
     for ( int j1 = 0; j1 < alignsx.isize( ); j1++ )
     for ( int j2 = j1 + 1; j2 < alignsx.isize( ); j2++ )
     {    if ( alignsx[j1].first != alignsx[j2].first ) continue;
          {    if ( Meet( alignsx[j1].second, alignsx[j2].second ) )
                    E.Join( j1, j2 );    }    }
     vec<int> reps, o;
     E.OrbitReps(reps);
     vec< pair< int, vec< pair<int,int> > > > alignsy;
     vec<double> mismatch_ratesy;
     START = -1, STOP = -1;
     for ( int j = 0; j < reps.isize( ); j++ )
     {    vec< pair<int,int> > M;
          E.Orbit( reps[j], o );
          Sort(o);
          for ( int l = 0; l < o.isize( ); l++ )
               M.append( alignsx[ o[l] ].second );
          UniqueSort(M);
          if ( BinMember( o, 0 ) ) START = alignsy.size( );
          if ( BinMember( o, 1 ) ) STOP = alignsy.size( );
          alignsy.push( alignsx[ o[0] ].first, M );    
          double mr = 1.0;
          for ( int l = 0; l < o.isize( ); l++ )
               mr = Min( mr, mismatch_ratesx[ o[l] ] );
          mismatch_ratesy.push_back(mr);    }
     alignsx = alignsy;
     mismatch_ratesx = mismatch_ratesy;

     // Screen alignments.
          
     pair< int, vec< pair<int,int> > > STARTA, STOPA;
     STARTA = alignsx[START], STOPA = alignsx[STOP];
     vec<Bool> to_remove( alignsx.size( ), False );
     SortSync( mismatch_ratesx, alignsx );
     vec<Bool> cov( target.size( ), False );
     for ( int i = 0; i < alignsx.isize( ); i++ )
     {    const vec< pair<int,int> >& M = alignsx[i].second;
          Bool more = False;
          for ( int j = 0; j < M.isize( ); j++ )
          {    int start = M.front( ).second, stop = M.back( ).second + L;
               for ( int l = start + K/2 - 1; l < stop - K/2; l++ )
               {    if ( !cov[l] ) more = True;
                    cov[l] = True;    }    }
          if ( alignsx[i] == STARTA || alignsx[i] == STOPA ) continue;
          if ( !more ) to_remove[i] = True;    }
     EraseIf( alignsx, to_remove );
     START = Position( alignsx, STARTA ), STOP = Position( alignsx, STOPA );

     // Print alignments.

     if ( verbosity >= 1 ) PRINT_TO( rout, alignsx.size( ) );
     if ( verbosity >= 2 )
     {    for ( int j = 0; j < alignsx.isize( ); j++ )
          {    rout << alignsx[j].first << ":";
               for ( int z = 0; z < alignsx[j].second.isize( ); z++ )
               {    int z2;
                    for ( z2 = z + 1; z2 < alignsx[j].second.isize( ); z2++ )
                    {    if ( alignsx[j].second[z2].first 
                              != alignsx[j].second[z2-1].first + 1 )
                         {    break;    }
                         if ( alignsx[j].second[z2].second 
                              != alignsx[j].second[z2-1].second + 1 )
                         {    break;    }    }
                    rout << " (" << alignsx[j].second[z].first << "-" 
                         << alignsx[j].second[z2-1].first + L << ","
                         << alignsx[j].second[z].second << "-" 
                         << alignsx[j].second[z2-1].second + L << ")";
                    z = z2 - 1;    }
               rout << endl;    }    }    }

void PickBestPath( const vec<basevector>& bpaths, const basevector& target,
     basevector& bpath, const int verbosity, ostream& rout, const String& run_dir,
     const String& data_dir )
{
     alignment al;
     int best_loc;
     vec<int> errs( bpaths.size( ) );
     for ( int j = 0; j < bpaths.isize( ); j++ )
     {    if ( bpaths[j].size( ) <= target.size( ) )
          {    SmithWatFree( bpaths[j], target, best_loc, al, True, True );
               errs[j] = ActualErrors( bpaths[j], target, align(al), 2, 3 );    }
          else
          {    SmithWatFree( target, bpaths[j], best_loc, al, True, True );
               errs[j] = ActualErrors(
                    target, bpaths[j], align(al), 2, 3 );    }    }
     vec<int> ids( errs.size( ), vec<int>::IDENTITY );
     SortSync( errs, ids );
     bpath = bpaths[ ids[0] ];
     if ( verbosity >= 2 )
     {    PRINT_TO( rout, ids[0] );
          rout << "\nbest bpath has " << errs[0] << " errors, "
               << "alignment to reference:\n";
          vecbasevector B(1);
          B[0] = bpath;
          #pragma omp critical
          {    AlignToRef( B, run_dir, data_dir, 
                    rout, "CleanLongReadPatches" );    }    }    }
