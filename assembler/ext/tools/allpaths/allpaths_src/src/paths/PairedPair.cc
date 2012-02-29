/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Set.h"
#include "VecOverlap.h"
#include "math/Functions.h"
#include "math/HoInterval.h"
#include "paths/FindClosures.h"
#include "paths/HyperKmerPath.h"
#include "paths/PairedPair.h"
#include "paths/SimpleWalk.h"

istring::display_mode istring::s_displayMode = istring::SPLIT;

int istring::Length( const vec<int>& L ) const
{    int len = 0;
     for ( int i = 0; i < isize( ); i++ )
          len += L[ (*this)[i] ];
     return len;    }

void GetOverlaps( const pp_read& r1, const pp_read& r2, vec<int>& offsets,
     const Bool require_first_subset_second )
{    offsets.clear( );
     for ( int o = -r2.isize( ) + 1; o <= r1.isize( ) - 1; o++ )
     {    if (require_first_subset_second)
               if ( o > 0 || o + r2.isize( ) < r1.isize( ) ) continue;
          int start1 = Max( 0, o ), stop1 = Min( r1.isize( ), r2.isize( ) + o );
          int j;
          for ( j = start1; j < stop1; j++ )
               if ( r1[j] != r2[j-o] ) break;
          if ( j == stop1 ) offsets.push_back(o);    }    }

void GetOverlaps( const pp_read& r1, const pp_read& r2, vec<int>& offsets,
     const vec<int>& L, int min_overlap )
{    offsets.clear( );
     for ( int o = -r2.isize( ) + 1; o <= r1.isize( ) - 1; o++ )
     {    int start1 = Max( 0, o ), stop1 = Min( r1.isize( ), r2.isize( ) + o );
          int j;
          int overlap = 0;
          for ( j = start1; j < stop1; j++ )
          {    if ( r1[j] != r2[j-o] ) break;
               overlap += L[ r1[j] ];    }
          if ( j == stop1 && overlap >= min_overlap ) offsets.push_back(o);    }    }

void JoinReads( const pp_read& r1, const pp_read& r2, const int offset, pp_read& r )
{    int left = Min( 0, offset ), right = Max( r1.isize( ), r2.isize( ) + offset );
     r.reserve( right - left );
     if ( offset >= 0 ) 
     {    r = r1;
          for ( int j = r1.isize( ) - offset; j < r2.isize( ); j++ )
               r.push_back( r2[j] );    }
     else
     {    r = r2;
          for ( int j = r2.isize( ) + offset; j < r1.isize( ); j++ )
               r.push_back( r1[j] );    }    }

int OverlapLength( const pp_read& r1, const pp_read& r2, const int offset, 
     const vec<int>& L )
{    int left = Max( 0, offset ), right = Min( r1.isize( ), r2.isize( ) + offset );
     int overlap = 0;
     for ( int j = left; j < right; j++ )
          overlap += L[ r1[j] ];
     return overlap;    }
          
int TransferMark1( const pp_read& r1, const pp_read& r2, const int offset, int m1 )
{    if ( offset < 0 ) m1 -= offset;
     return m1;    }

int TransferMark2( const pp_read& r1, const pp_read& r2, const int offset, int m2 )
{    if ( offset > 0 ) m2 += offset;
     return m2;    }

Bool ReadyToClose( const pp_pair& p, const double dmult )
{    return p.Gap( ) < -dmult * p.Dev( );    }

// FindJoins: Given two pp_reads r1, r2 whose positions are given by
// p1 +/- d1, p +/- d2, respectively, try to join them, yielding new pp_reads
// ri having positions pi +/- di.

void FindJoins( 
     /* inputs: */  const pp_read& r1, const double p1, const double d1, 
                    const pp_read& r2, const double p2, const double d2, 
                    const double dmult, const vec<int>& L, 
     /* outputs: */ vec<pp_read>& r, vec<double>& p, vec<double>& d, 
                    vec<int>& offsets )
{    r.clear( ), p.clear( ), d.clear( ), offsets.clear( );
     static vec<int> offsets0;
     GetOverlaps( r1, r2, offsets0 );
     for ( int i = 0; i < offsets0.isize( ); i++ )
     {    int o = offsets0[i];
          double p1x(p1), p2x(p2);
          if ( o > 0 )
          {    for ( int j = 0; j < o; j++ )
                    p2x -= L[ r1[j] ];    }
          if ( o < 0 )
          {    for ( int j = 0; j < -o; j++ )
                    p1x -= L[ r2[j] ];    }
          double px, dx;
          if ( CombineMeasurements( p1x, p2x, d1, d2, dmult, px, dx ) )
          {    static pp_read rx;
               JoinReads( r1, r2, o, rx );
               r.push_back(rx), p.push_back(px), d.push_back(dx);
               offsets.push_back(o);    }    }    }

// FindJoins- Warning!  If wt1 and wt2 are positive, their values are IGNORED!

void FindJoins( const pp_mpair& p, const pp_pair& x, const vec<int>& L, 
     const double dmult, vec<pp_mpair>& pnew, int wt1, int wt2, Bool must_extend )
{    pnew.clear( );
     vec<int> o1, o2;
     GetOverlaps( p.Left( ), x.Left( ), o1 );
     if ( o1.empty( ) ) return;
     GetOverlaps( p.Right( ), x.Right( ), o2 );
     if ( o2.empty( ) ) return;
     if (must_extend)
     {    if ( o1.size( ) == 1 && o1[0] >= 0 
               && o1[0] + x.LeftSize( ) <= p.LeftSize( )
               && o2.size( ) == 1 && o2[0] >= 0 
               && o2[0] + x.RightSize( ) <= p.RightSize( ) )
          {    return;    }    }
     for ( int i1 = 0; i1 < o1.isize( ); i1++ )
     {    pp_read y1, y2;
          JoinReads( p.Left( ), x.Left( ), o1[i1], y1 );
          int m1 = TransferMark1( p.Left( ), x.Left( ), o1[i1], p.LeftMark( ) );
          for ( int i2 = 0; i2 < o2.isize( ); i2++ )
          {    JoinReads( p.Right( ), x.Right( ), o2[i2], y2 );
               int m2 = TransferMark1( 
                    p.Right( ), x.Right( ), o2[i2], p.RightMark( ) );
               double g1 = p.Gap( ), g2 = x.Gap( );
               for ( int j = p.LeftSize( ) - o1[i1]; j < x.LeftSize( ); j++ )
                    g1 -= L[ x.Left(j) ];
               for ( int j = 0; j < -o2[i2]; j++ )
                    g1 -= L[ x.Right(j) ];
               for ( int j = o1[i1] + x.LeftSize( ); j < p.LeftSize( ); j++ )
                    g2 -= L[ p.Left(j) ];
               for ( int j = 0; j < o2[i2]; j++ )
                    g2 -= L[ p.Right(j) ];
               if ( wt1 > 0 && wt2 > 0 )
               {    double d1 = p.Dev( ), d2 = x.Dev( );
                    double g, d;
                    if ( !CombineMeasurements( g1, g2, d1, d2, dmult, g, d ) )
                         continue;
                    pnew.push_back( pp_mpair( y1, y2, m1, m2, g, d ) );    }
               else
               {    if ( Abs( g1 - g2 ) > dmult * ( p.Dev( ) + x.Dev( ) ) ) continue;
                    pnew.push_back( pp_mpair( y1, y2, m1, m2, 
                         g1, p.Dev( ) ) );    }    }    }
     UniqueSort(pnew);    }

void FindExtensions( const pp_mpair& p, const vec<pp_pair>& pairs, 
     const vec<int>& L, const double dmult, vec<pp_mpair>& pnew )
{    pnew.clear( );
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    const pp_pair& x = pairs[i];
          static vec<pp_mpair> pnewi;
          FindJoins( p, x, L, dmult, pnewi );
          pnew.append(pnewi);    }
     UniqueSort(pnew);    }

void GetClosures( const pp_mpair& p, const vec<int>& L, const double dmult, 
     vec<pp_closure>& closures, const Bool trim, const Bool verbose )
{    closures.clear( );
     vec<int> offsets;
     GetOverlaps( p.Left( ), p.Right( ), offsets );
     for ( int j = 0; j < offsets.isize( ); j++ )
     {    int o = offsets[j];
          int overlap = 0;
          for ( int k = Max( 0, o ); k < p.LeftSize( ); k++ )
               overlap += L[ p.Left(k) ];
          for ( int k = 0; k < -o; k++ )
               overlap += L[ p.Right(k) ];
          double dm = dmult;
          for ( int pass = 1; pass <= 2; pass++ )
          {    if ( pass == 2 )
               {    if ( closures.nonempty( ) ) break;
                    dm += 1.0;    }
               double offby = Abs( double(overlap) + p.Gap( ) );
               if ( offby <= dm * p.Dev( ) )
               {    if (verbose) 
                    {    cout << "see closure with dev = ";
                         if ( p.Dev( ) > 0 )
                              cout << setprecision(3) << offby / p.Dev( );
                         else if ( p.Dev( ) == 0.0 && offby == 0.0 ) cout << 0;
                         else if ( p.Dev( ) == 0.0 && offby > 0 ) cout << "infinity";
                         cout << "\n";    }
                    pp_read r;
                    JoinReads( p.Left( ), p.Right( ), o, r );
                    int m1 = TransferMark1( p.Left(), p.Right(), o, p.LeftMark() );
                    int m2 = TransferMark2( p.Left(), p.Right(), o, p.RightMark() );
                    pp_closure c;
                    if ( !trim ) c = r;
                    else
                    {    c.clear( );
                         for ( int k = m1; k < m2; k++ )
                              c.push_back( r[k] );    }
                    closures.push_back(c);    }    }    }
     UniqueSort(closures);    }

// ImpliedGap.  Suppose given a read x and alignments (x,pL) and (x,pR) with
// respective offsets lo and ro.  Return the gap that is implied between pL and pR.
// The vector L is to contain the lengths.

int ImpliedGap( 
     const pp_read& x, const pp_pair& p, const vec<int>& L, int lo, int ro )
{    int gap = 0;
     if ( lo + p.LeftSize( ) <= ro )
     {    for ( int j = lo + p.LeftSize( ); j < ro; j++ )
               gap += L[ x[j] ];
          return gap;    }
     for ( int j = Max( 0, ro ); j < Min( x.isize( ), lo + p.LeftSize( ) ); j++ )
          gap -= L[ x[j] ];
     for ( int j = 0; j < -ro; j++ )
          gap -= L[ p.Right(j) ];
     for ( int j = x.isize( ); j < lo + p.LeftSize( ); j++ )
          gap -= L[ p.Left(j-lo) ];
     return gap;    }

int ImpliedUnitGap( const pp_read& x, const pp_pair& p, int lo, int ro )
{    int gap = 0;
     if ( lo + p.LeftSize( ) <= ro )
     {    for ( int j = lo + p.LeftSize( ); j < ro; j++ )
               gap++;
          return gap;    }
     for ( int j = Max( 0, ro ); j < Min( x.isize( ), lo + p.LeftSize( ) ); j++ )
          gap--;
     for ( int j = 0; j < -ro; j++ )
          gap--;
     for ( int j = x.isize( ); j < lo + p.LeftSize( ); j++ )
          gap--;
     return gap;    }

// ValidOffset: given an offset o between reads (x,y), return True unless the offset
// forces them to overlap but they don't agree.

Bool ValidOffset( const pp_read& x, const pp_read& y, int o )
{    for ( int j = Max( 0, -o ); j < Min( y.isize( ), x.isize( ) - o ); j++ )
          if ( x[j+o] != y[j] ) return False;
     return True;    }

// ValidPairOverlap: given two pairs p and q, overlapping on the left with offset
// llo and on the right with offset rro, are the overlaps consistent with the 
// gaps for p and q and the possible overlaps (pL,pR) as defined by poffsets and
// the possible overlaps (qL,qR) as defined by qoffsets?

Bool ValidPairOverlap( const pp_pair& p, const pp_pair& q, int llo, int rro,
     const vec<int>& poffsets, const vec<int>& qoffsets, const vec<int> L,
     double dmult )
{
     // First suppose that pL overlaps pR.

     double pqdev = p.Dev( ) + q.Dev( );
     for ( int i = 0; i < poffsets.isize( ); i++ )
     {    int po = poffsets[i];
          int pgap = ImpliedGap( p.Left( ), p, L, 0, po );
          if ( Abs( p.Gap( ) - pgap ) > dmult * pqdev ) continue;
          int qo = po + rro - llo;
          int qgap = 0, qgapu = 0;
          if ( qo <= -q.RightSize( ) ) continue;
          if ( qo <= q.LeftSize( ) )
          {    qgap = ImpliedGap( q.Left( ), q, L, 0, qo );
               qgapu = ImpliedUnitGap( q.Left( ), q, 0, qo );    }
          else
          {    for ( int j = q.LeftSize( ); j < qo; j++ )
               {    if ( j + llo < p.LeftSize( ) ) qgap += L[ p.Left(j+llo) ];
                    else qgap += L[ p.Right(j+llo-po) ];
                    qgapu++;    }    }
          if ( Abs( q.Gap( ) - qgap ) > dmult * pqdev ) continue;
          if ( !ValidOffset( q.Left( ), q.Right( ), q.LeftSize( ) + qgapu ) )
               continue;
          if ( !ValidOffset( p.Right( ), q.Left( ), llo - po ) ) continue;
          if ( !ValidOffset( p.Left( ), q.Right( ), llo + q.LeftSize( ) + qgapu ) )
               continue;
          return True;    }

     // Now suppose that qL overlaps qR.

     for ( int i = 0; i < qoffsets.isize( ); i++ )
     {    int qo = qoffsets[i];
          int qgap = ImpliedGap( q.Left( ), q, L, 0, qo );
          if ( Abs( q.Gap( ) - qgap ) > dmult * pqdev ) continue;
          int po = qo + llo - rro;
          int pgap = 0, pgapu = 0;
          if ( po <= -p.RightSize( ) ) continue;
          if ( po <= p.LeftSize( ) ) 
          {    pgap = ImpliedGap( p.Left( ), p, L, 0, po );
               pgapu = ImpliedUnitGap( p.Left( ), p, 0, po );    }
          else
          {    for ( int j = p.LeftSize( ); j < po; j++ )
               {    if ( j < llo + q.LeftSize( ) ) pgap += L[ q.Left(j-llo) ];
                    else pgap += L[ q.Right(j-llo-qo) ];
                    pgapu++;    }    }
          if ( Abs( p.Gap( ) - pgap ) > dmult * pqdev ) continue;
          if ( !ValidOffset( p.Left( ), p.Right( ), p.LeftSize( ) + pgapu ) )
               continue;
          if ( !ValidOffset( p.Left( ), q.Right( ), llo + qo ) ) continue;
          if ( !ValidOffset( p.Right( ), q.Left( ), llo - p.LeftSize( ) - pgapu ) )
               continue;
          return True;    }

     // Now we can assume that there is no overlap on either side.

     int gapdiff = 0; // will be pgap - qgap
     for ( int i = p.LeftSize( ); i < llo + q.LeftSize( ); i++ )
          gapdiff += L[ q.Left(i-llo) ];
     for ( int i = llo + q.LeftSize( ); i < p.LeftSize( ); i++ )
          gapdiff -= L[ p.Left(i) ];
     for ( int i = 0; i < -rro; i++ )
          gapdiff += L[ q.Right(i) ];
     for ( int i = 0; i < rro; i++ )
          gapdiff -= L[ p.Right(i) ];
     double epsilon = 0.000001;
     return IntervalOverlap(
          Max( 0.0, p.Gap( ) - dmult * p.Dev( ) - epsilon ), 
          p.Gap( ) + dmult * p.Dev( ),
          Max( 0.0, q.Gap( ) - dmult * q.Dev( ) - epsilon ) + gapdiff, 
          q.Gap( ) + dmult * q.Dev( ) + epsilon + gapdiff ) > 0;    }

// ValidateClosures.  For each closure c, first find all left and right extensions 
// by reads, gluing on to yield closures c'.  Then for each c', find placements of 
// pairs on it, with both sides of the pair hitting c'.  Then determine if one can 
// cover the gap in the original pair by iteratively layering on the placed pairs.
//
// This version will not see pairs that hit c' in the following way:
//      ------left extension-------...........c.............
//         =====pL====        ===pR======
// where pR reaches off the left end of c.
//
// Added feature:
//
// Pretest 1.  We first attempt to validate a closure c of a pair p by looking 
// for another pair q so that we have the following picture:
//
//              ........................c........................
//              =========pL========             =======pR========
//          ----------------qL--------------------                   -----qR------
//
// or the reverse of it:
//
//              ........................c........................
//              =========pL========             =======pR========
//  ----qL---                    --------------qR--------------------

void ValidateClosures( const pp_pair& p, vec<pp_closure>& closures, 
     const vec<pp_pair>& pairs, const vec_overlap<int>& over,
     const vec< vec<int> > reads, const vec<int>& L, const double dmult,
     const int verbosity )
{    static vec<Bool> to_delete;
     to_delete.resize_and_set( closures.size( ), False );
     vec<Bool> Validated( closures.size( ), False );
     for ( int z = 0; z < closures.isize( ); z++ )
     {    const pp_closure& c = closures[z];
          static vec< pair<int,int> > overlaps;
          over.GetOverlaps( c, overlaps );

          // Pretest 1.

          for ( int v = 0; v < overlaps.isize( ); v++ )
          {    int id = overlaps[v].first, o = overlaps[v].second;
               const pp_pair& q = pairs[id/2];
               int c_extend = 0;
               if ( id % 2 == 0 )
               {    if ( o + q.LeftSize( ) <= c.isize( ) - p.RightSize( ) ) continue;
                    if ( o >= p.LeftSize( ) ) continue;
                    if ( o + q.LeftSize( ) < c.isize( ) )
                    {    for ( int j = o + q.LeftSize( ); j < c.isize( ); j++ )
                              c_extend += L[ c[j] ];
                         if ( q.Gap( ) - c_extend < -dmult * q.Dev( ) ) 
                              continue;    }
                    // cout << "\n1validating closure " << c << "\nof " << p
                    //      << "\nusing " << q << "\n";
                    Validated[z] = True;
                    break;    }
               else
               {    if ( o >= p.LeftSize( ) ) continue;
                    if ( o + q.RightSize( ) <= c.isize( ) - p.RightSize( ) ) 
                         continue;
                    if ( o > 0 )
                    {    for ( int j = 0; j < o; j++ )
                              c_extend += L[ c[j] ];
                         if ( q.Gap( ) - c_extend < -dmult * q.Dev( ) ) 
                              continue;    }
                    // cout << "\n2validating closure " << c << "\nof " << p
                    //      << "\nusing " << q << "\n";
                    Validated[z] = True;
                    break;    }    }
          if ( Validated[z] ) continue;

          // Continue with main test.

          static vec<pp_read> left_exts, right_exts;
          left_exts.clear( ), right_exts.clear( );
          static pp_read empty_read;
          left_exts.push_back(empty_read), right_exts.push_back(empty_read);
          for ( int v = 0; v < overlaps.isize( ); v++ )
          {    int id = overlaps[v].first, o = overlaps[v].second;
               static pp_read e;
               const vec<int>& r = reads[id];
               if ( id % 2 == 0 && o < 0 )
               {    e.SetToRangeOf( r, 0, -o );
                    left_exts.push_back(e);    }
               if ( id % 2 == 1 && o + r.isize( ) > c.isize( ) )
               {    e.SetToRangeOf( r, c.isize( ) - o, r.isize( ) );
                    right_exts.push_back(e);    }    }
          UniqueSort(left_exts), UniqueSort(right_exts);
          int np = pairs.size( );
          vec< vec<int> > left_placements(np), right_placements(np);
          static vec<int> offsets;
          for ( int jl = 0; jl < left_exts.isize( ); jl++ )
          {    const pp_read& e = left_exts[jl];
               for ( int u = 0; u < pairs.isize( ); u++ )
               {    const pp_pair& P = pairs[u];
                    GetOverlaps( e, P.Left( ), offsets );
                    for ( int m1 = 0; m1 < offsets.isize( ); m1++ )
                    {    int stop1 = offsets[m1] + P.LeftSize( );
                         if ( stop1 > e.isize( ) ) continue;
                         int dist = 0;
                         for ( int i = stop1; i < e.isize( ); i++ )
                              dist += L[ e[i] ];
                         left_placements[u].push_back(dist);    }    }    }
          for ( int jl = 0; jl < right_exts.isize( ); jl++ )
          {    const pp_read& e = right_exts[jl];
               for ( int u = 0; u < np; u++ )
               {    const pp_pair& P = pairs[u];
                    GetOverlaps( e, P.Right( ), offsets );
                    for ( int m2 = 0; m2 < offsets.isize( ); m2++ )
                    {    int start2 = offsets[m2];
                         if ( start2 < 0 ) continue;
                         int dist = 0;
                         for ( int i = 0; i < start2; i++ )
                              dist += L[ e[i] ];
                         right_placements[u].push_back(dist);    }    }    }
          for ( int u = 0; u < pairs.isize( ); u++ )
          {    UniqueSort( left_placements[u] );
               UniqueSort( right_placements[u] );    }
          vec< vec<int> > middle_left_placements(np), middle_right_placements(np);
          for ( int v = 0; v < overlaps.isize( ); v++ )
          {    int id = overlaps[v].first, o = overlaps[v].second;
               static pp_read e;
               const vec<int>& r = reads[id];
               if ( id % 2 == 0 )
                    middle_left_placements[id/2].push_back(o);
               else middle_right_placements[id/2].push_back(o);    }
          static vec<ho_interval> covers;
          covers.clear( );
          for ( int u = 0; u < pairs.isize( ); u++ )
          {    const pp_pair& P = pairs[u];
               for ( int j1 = 0; j1 < middle_left_placements[u].isize( ); j1++ )
               {    int o1 = middle_left_placements[u][j1];
                    for ( int j2 = 0; j2 < middle_right_placements[u].isize(); j2++ )
                    {    int o2 = middle_right_placements[u][j2];
                         int start1 = o1, stop1 = o1 + P.LeftSize( );
                         start1 = Max( 0, start1 );
                         stop1 = Min( c.isize( ), stop1 );
                         int start2 = o2, stop2 = o2 + P.RightSize( );
                         start2 = Max( 0, start2 );
                         stop2 = Min( c.isize( ), stop2 );
                         int X1 = o1 + P.LeftSize( ), X2 = o2;
                         int gap = 0;
                         if ( X1 <= X2 )
                         {    for ( int s = X1; s < X2; s++ )
                                   gap += L[ c[s] ];    }
                         else
                         {    if ( o2 - o1 < 0 ) continue;
                              if ( o1 + P.LeftSize( ) > o2 + P.RightSize( ) )
                                   continue;
                              for ( int s = o2 - o1; s < P.LeftSize( ); s++ )
                                   gap -= L[ P.Left(s) ];    }
                         if ( Abs( gap - P.Gap( ) ) > dmult * P.Dev( ) ) continue;
                         covers.push_back( ho_interval( start1, stop1 ) );
                         covers.push_back( ho_interval( start2, stop2 ) );    }
                    for ( int j2 = 0; j2 < right_placements[u].isize(); j2++ )
                    {    int d2 = right_placements[u][j2];
                         int start1 = o1, stop1 = o1 + P.LeftSize( );
                         start1 = Max( 0, start1 );
                         stop1 = Min( c.isize( ), stop1 );
                         int X1 = o1 + P.LeftSize( ), X2 = c.size( );
                         if ( X1 <= X2 )
                         {    int gap = d2;
                              for ( int s = X1; s < X2; s++ )
                                   gap += L[ c[s] ];
                              if ( Abs( gap - P.Gap( ) ) > dmult * P.Dev( ) )
                                   continue;
                              covers.push_back( ho_interval( 
                                   start1, stop1 ) );    }    }    }
               for ( int j1 = 0; j1 < left_placements[u].isize( ); j1++ )
               {    int d1 = left_placements[u][j1];
                    for ( int j2 = 0; j2 < middle_right_placements[u].isize(); j2++ )
                    {    int o2 = middle_right_placements[u][j2];
                         int start2 = o2, stop2 = o2 + P.RightSize( );
                         start2 = Max( 0, start2 );
                         stop2 = Min( c.isize( ), stop2 );
                         int X1 = 0, X2 = o2;
                         if ( X1 <= X2 )
                         {    int gap = d1;
                              for ( int s = X1; s < X2; s++ )
                                   gap += L[ c[s] ];
                              if ( Abs( gap - P.Gap( ) ) > dmult * P.Dev( ) )
                                   continue;
                              covers.push_back( ho_interval( 
                                   start2, stop2 ) );    }    }    }    }

          Sort( covers );
          if ( covers.empty() || covers.front().Start() != 0 ) 
            to_delete[z] = True;
          else {
            int maxStop = covers.front().Stop();
            for ( int cIdx = 1; cIdx < covers.isize(); ++cIdx )
              if ( covers[cIdx].Start() >= maxStop ) {
                to_delete[z] = True;
                break;
              } else if ( covers[cIdx].Stop() > maxStop )
                maxStop = covers[cIdx].Stop();
            if ( maxStop < c.isize() )
              to_delete[z] = True;
          }
     }          

     EraseIf( closures, to_delete );
     EraseIf( Validated, to_delete );

     // Phase 2: make sure the the same left and right extension work for the
     // given set of closures.

     to_delete.resize_and_set( closures.size( ), False );
     for ( int z = 0; z < closures.isize( ); z++ )
     {    if ( Validated[z] ) continue;
          const pp_closure& c = closures[z];
          if ( verbosity > 2 ) PRINT( c );
          static vec<pp_read> left_exts, right_exts;
          left_exts.clear( ), right_exts.clear( );
          static pp_read empty_read;
          left_exts.push_back(empty_read), right_exts.push_back(empty_read);
          for ( int u = 0; u < pairs.isize( ); u++ )
          {    for ( int pass = 1; pass <= 2; pass++ )
               {    const pp_read& r 
                         = ( pass == 1 ? pairs[u].Left( ) : pairs[u].Right( ) );
                    static vec<int> offsets;
                    GetOverlaps( c, r, offsets );
                    for ( int v = 0; v < offsets.isize( ); v++ )
                    {    int o = offsets[v];
                         static pp_read e;
                         if ( o < 0 )
                         {    e.SetToRangeOf( r, 0, -o );
                              left_exts.push_back(e);    }
                         else if ( o + r.isize( ) > c.isize( ) )
                         {    e.SetToRangeOf( r, c.isize( ) - o, r.isize( ) );
                              right_exts.push_back(e);    }    }    }    }
          UniqueSort(left_exts), UniqueSort(right_exts);
          Bool validated = False;
          for ( int jl = 0; jl < left_exts.isize( ); jl++ )
          {    if (validated) break;
               for ( int jr = 0; jr < right_exts.isize( ); jr++ )
               {    if (validated) break;
                    if ( verbosity > 3 ) {
                      PRINT2( left_exts[jl], right_exts[jr] );
                      for ( int z = 0; z < left_exts[jl].isize(); ++z )
                        cout << "-";
                      for ( int z = 0; z < p.LeftSize()-1; ++z )
                        cout << "=";
                      cout << ">";
                      for (int z = p.LeftSize() + p.RightSize(); z < c.isize(); ++z)
                        cout << "=";
                      cout << "<";
                      for ( int z = 0; z < p.RightSize()-1; ++z )
                        cout << "=";
                      for ( int z = 0; z < right_exts[jr].isize(); ++z )
                        cout << "-";
                      cout << "\t\t" << p.Left() << "   " << p.Right() << endl;
                    }
                    static pp_closure cp;
                    cp = left_exts[jl];
                    cp.append(c);
                    cp.append( right_exts[jr] );
                    int left = left_exts[jl].size( ) + p.LeftSize( ); 
                    int right 
                         = cp.isize( ) - right_exts[jr].isize( ) - p.RightSize( );
                    static vec< pair<ho_interval,ho_interval> > covers;
                    covers.clear( );
                    for ( int u = 0; u < pairs.isize( ); u++ )
                    {    const pp_pair& P = pairs[u];
                         static vec<int> offsets1, offsets2;
                         GetOverlaps( P.Left( ), cp, offsets1 );
                         GetOverlaps( P.Right( ), cp, offsets2 );
                         for ( int m1 = 0; m1 < offsets1.isize( ); m1++ )
                         {    int o1 = -offsets1[m1];
                              int start1 = o1, stop1 = o1 + P.LeftSize( );
                              start1 = Max( 0, start1 );
                              stop1 = Min( cp.isize( ), stop1 );
                              for ( int m2 = 0; m2 < offsets2.isize( ); m2++ )
                              {    int o2 = -offsets2[m2];
                                   int start2 = o2;
                                   int stop2 = o2 + P.RightSize( );
                                   start2 = Max( 0, start2 );
                                   stop2 = Min( cp.isize( ), stop2 );
                                   int X1 = o1 + P.LeftSize( ), X2 = o2;

                                   int gap = 0;
                                   if ( X1 <= X2 )
                                   {    for ( int s = X1; s < X2; s++ )
                                             gap += L[ cp[s] ];    }
                                   else
                                   {    if ( o2 - o1 < 0 ) continue;
                                        if (o1 + P.LeftSize( ) > o2 + P.RightSize( ))
                                             continue;
                                        for (int s = o2 - o1; s < P.LeftSize( ); s++)
                                             gap -= L[ P.Left(s) ];    }
                                   if ( Abs( gap - P.Gap( ) ) > dmult * P.Dev( ) ) 
                                        continue;
                                   if ( verbosity > 3 ) {
                                     if ( gap > 0 || P.Left() != P.Right() ) {
                                       for ( int z = 0; z < o1; ++z )
                                         cout << " ";
                                       for ( int z = start1; z < stop1-1; ++z )
                                         cout << "-";
                                       cout << ">";
                                       for ( int z = stop1; z < start2; ++z )
                                         cout << " ";
                                       cout << "<";
                                       for ( int z = start2+1; z < stop2; ++z )
                                         cout << "-";
                                       cout << "\t\t" << P.Left() 
                                            << "   " << P.Right() << endl;
                                     } else {
                                       cout << "<";
                                       for ( int z = start1+1; z < stop1-1; ++z )
                                         cout << "-";
                                       cout << ">";
                                       cout << "\t\t" << P.Left() << endl;
                                     }
                                   }
                                   covers.push_back( make_pair( 
                                        ho_interval( start1, stop1 ),
                                        ho_interval( 
                                             start2, stop2 ) ) );    }    }    }
                    while( left <= right )
                    {    Bool progress = False;
                         for ( int u = 0; u < covers.isize( ); u++ )
                         {    const ho_interval &h1 = covers[u].first; 
                              const ho_interval &h2 = covers[u].second;
                              if ( h1.Start( ) >= left ) continue;
                              if ( h2.Stop( ) <= right ) continue;
                              if ( h1.Stop( ) > left || h2.Start( ) < right )
                              {    progress = True;
                                   left = Max( left, h1.Stop( ) );
                                   right = Min( right, h2.Start( ) );    }    }
                         if ( !progress ) break;    }
                    if ( left > right ) {
                      validated = True; 
                      if ( verbosity > 2 )
                        cout << "validated" << endl;
                    }    }    }
          if ( !validated ) to_delete[z] = True;    }
     EraseIf( closures, to_delete );    }

void ClosePairedPairsEasy( int min_side, const HyperKmerPath& h, 
     const vec<pp_pair>& pairs, const vec<Bool>& pairs_to_close, const vec<int>& L, 
     const double dmult, vec< vec<pp_closure> >& closures, vec< vec<double> >& devs,
     int max_ext, vec<Bool>& fail, int max_processed, int max_unprocessed, 
     int verbosity, int simple_walk_verbosity, int max_opens, int max_nodes,
     Bool create_closures_if_fail, const int MAX_CLOSURES )
{    
     if ( verbosity >= 1 ) cout << "\nClosePairedPairsEasy:\n";
     double eclock = 0.0;
     if ( verbosity >= 1 ) eclock = WallClockTime( );
     closures.clear( ), devs.clear( );
     closures.resize( pairs.size( ) );
     devs.resize( pairs.size( ) );
     fail.resize_and_set( pairs.size( ), True );
     static vec<int> leftsize, rightsize;
     leftsize.resize_and_set( pairs.size( ), 0 );
     rightsize.resize_and_set( pairs.size( ), 0 );
     vec< vec<int> > readsp;
     for ( int u = 0; u < pairs.isize( ); u++ )
          readsp.push_back( pairs[u].Left( ), pairs[u].Right( ) );
     vec_overlap<int> over(readsp);
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    const pp_pair& p = pairs[i];
          for ( int j = 0; j < p.LeftSize( ); j++ )
               leftsize[i] += L[ p.Left(j) ];
          for ( int j = 0; j < p.RightSize( ); j++ )
               rightsize[i] += L[ p.Right(j) ];    }
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    if ( !pairs_to_close[i] ) continue;
          const pp_pair& p = pairs[i];
          double pdev_max = dmult * p.Dev( );
          if ( leftsize[i] < min_side && rightsize[i] < min_side ) continue;

          // Not sure if we want to do the cases where the pair closes itself, but
          // for now, we skip them:

          if ( p.Gap( ) < -dmult * p.Dev( ) ) continue;

          if ( verbosity >= 1 && leftsize[i] >= min_side )
          {    cout << "\n--------------------------------------------------"
                    << "----------------------------------\n";
               cout << "\nTo close " << p << "\n\n";    }
          if ( leftsize[i] >= min_side )
          {    fail[i] = False;
               static vec<int> use;
               use.clear( );
               for ( int j = 0; j < pairs.isize( ); j++ )
               {    if ( j == i ) continue;
                    const pp_pair& q = pairs[j];
                    static vec<int> lloffsets, rroffsets, lroffsets;
                    static vec<int> poffsets, qoffsets;
                    GetOverlaps( p.Left( ), q.Left( ), lloffsets );
                    if ( lloffsets.empty( ) ) continue;
                    GetOverlaps( p.Right( ), q.Right( ), rroffsets );
                    GetOverlaps( p.Left( ), q.Right( ), lroffsets );
                    GetOverlaps( p.Left( ), p.Right( ), poffsets );
                    GetOverlaps( q.Left( ), q.Right( ), qoffsets );
               
                    // Go through each possible overlap of pL and qL.

                    for ( int lu = 0; lu < lloffsets.isize( ); lu++ )
                    {    int llo = lloffsets[lu];

                         // First suppose that pR overlaps qR.

                         Bool hit = False;
                         for ( int l = 0; l < rroffsets.isize( ); l++ )
                         {    int rro = rroffsets[l];
                              if ( llo + q.LeftSize( ) > p.LeftSize( ) || rro < 0 )
                              {    if ( ValidPairOverlap( p, q, llo, rro, poffsets, 
                                        qoffsets, L, dmult ) )
                                   {    use.push_back(j);
                                        hit = True;
                                        break;    }    }    }
                         if (hit) break;

                         // From now on we can assume that pR does not overlap qR.
                         // Make sure that p's gap is big enough.

                         double g = p.Gap( ) + pdev_max - rightsize[j];
                         double pqdev = p.Dev( ) + q.Dev( );
                         for ( int k = p.LeftSize( ); k < llo + q.LeftSize( ); k++ )
                              g -= L[ q.Left( k - llo ) ];
                         for ( int k = llo + q.LeftSize( ); k < p.LeftSize( ); k++ )
                              g += L[ p.Left(k) ];
                         g -= ( q.Gap( ) - dmult * pqdev );
                         if ( g < 0 ) continue;

                         // Suppose that qR does not overlap pL.  We have to make
                         // sure that q's gap is big enough.

                         int poverhang = 0;
                         for ( int k = llo + q.LeftSize( ); k < p.LeftSize( ); k++ )
                              poverhang += L[ p.Left(k) ];
                         if ( poverhang <= q.Gap( ) + dmult * pqdev )
                         {    use.push_back(j);
                              break;    }

                         // Finally suppose that qR overlaps pL, extending it to the
                         // right.

                         for ( int lru = 0; lru < lroffsets.isize( ); lru++ )
                         {    int lro = lroffsets[lru];
                              if ( lro + q.RightSize( ) <= p.LeftSize( ) ) continue;
                              int qgap = ImpliedGap( p.Left( ), q, L, llo, lro );
                              if ( Abs( qgap - q.Gap( ) ) <= dmult * pqdev )
                              {    use.push_back(j);
                                   hit = True;
                                   break;    }    }
                         if (hit) break;    }    }

               if ( verbosity >= 2 )
               {    for ( int j = 0; j < use.isize( ); j++ )
                         cout << "[" << j+1 << "] " << pairs[ use[j] ] << "\n";    }

               static vec<pp_pair> pairsx;
               pairsx.clear( );
               pairsx.push_back( pairs[i] );
               for ( int j = 0; j < use.isize( ); j++ )
                    pairsx.push_back( pairs[ use[j] ] );

               /*
               static vec<Bool> pairs_to_closex;
               pairs_to_closex.clear( );
               pairs_to_closex.push_back(True);
               for ( int j = 0; j < use.isize( ); j++ )
                    pairs_to_closex.push_back(False);
               // BinaryOverwrite( "xxx.pairs", pairsx ); // BBBBBBBBBBBBBBBBBBBBBBB
               // BinaryOverwrite( "xxx.L", L ); // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               const int FINDCLOSURES_MAX_PSEUDO_CLOSURES = 1000000;
               const int FINDCLOSURES_MAX_CLOSURES = 1000000;
               vec< vec<pp_closure> > ppclosuresx;
               vec< vec<double> > devsx;
               vec<Bool> failx;
               unsigned int max_pseudo_closures = FINDCLOSURES_MAX_PSEUDO_CLOSURES;
               unsigned int max_closures = FINDCLOSURES_MAX_CLOSURES;
               double fclock = WallClockTime( );
               FindClosures( pairsx, pairs_to_closex, h, dmult, 
                    ppclosuresx, devsx, failx, max_pseudo_closures, max_closures );
               cout << TimeSince(fclock) << " used in FindClosures" << endl;
               cout << "\nfound " << ppclosuresx[0].size( ) << " closures\n";
               // Sort( ppclosuresx[0] ); // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               // for ( int j = 0; j < ppclosuresx[0].isize( ); j++ ) // BBBBBBBBBBB
               //      cout << "[" << j << "] " << ppclosuresx[0][j] << "\n"; // BBB
               // static int count(0); // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               // count++; // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               // if ( count == 2 ) exit(0); // BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
               */

               static vec<pp_read> reads;
               reads.clear( );
               for ( int j = 0; j < pairsx.isize( ); j++ )
                    reads.push_back( pairsx[j].Left( ), pairsx[j].Right( ) );
               UniqueSort(reads);
               static vec<pp_closure> ppclosuresx;
               Bool wfail;
               double sclock = WallClockTime( );
               SimpleWalkRight( pairs[i], reads, L, dmult, ppclosuresx, max_opens, 
                    create_closures_if_fail, wfail, simple_walk_verbosity,
                    False, max_nodes );
               if ( verbosity >= 1 )
               {    cout << TimeSince(sclock) << " used in SimpleWalkRight" << endl;
                    cout << "\nfound " << ppclosuresx.size( ) << " closures\n";    }
               if ( MAX_CLOSURES > 0 && ppclosuresx.isize( ) > MAX_CLOSURES )
               {    if ( verbosity >= 1 ) cout << "Too many closures to keep.\n";
                    fail[i] = True;
                    continue;    }
               if (wfail) fail[i] = True;

               // Validate the closures.

               double vclock = WallClockTime( );
               ValidateClosures( p, ppclosuresx, pairs, over, readsp, L, dmult, 
                    verbosity );
               if ( verbosity >= 1 )
               {    cout << TimeSince(vclock) << " used validating closures"
                         << endl;    }
               int valids = ppclosuresx.size( );
               for ( int j = 0; j < ppclosuresx.isize( ); j++ )
               {    const pp_closure& c = ppclosuresx[j];
                    double len = -pairs[i].Gap( );
                    for ( int u = 0; u < c.isize( ); u++ )
                         len += L[ c[u] ];
                    for ( int u = 0; u < pairs[i].LeftSize( ); u++ )
                         len -= L[ pairs[i].Left(u) ];
                    for ( int u = 0; u < pairs[i].RightSize( ); u++ )
                         len -= L[ pairs[i].Right(u) ];
                    double d = len / pairs[i].Dev( );
                    closures[i].push_back(c);    
                    devs[i].push_back(d);    }

               /*
               int valids = 0;
               for ( int j = 0; j < ppclosuresx.isize( ); j++ )
               {    const pp_closure& c = ppclosuresx[j];
                    double len = -pairs[i].Gap( );
                    for ( int u = 0; u < c.isize( ); u++ )
                         len += L[ c[u] ];
                    for ( int u = 0; u < pairs[i].LeftSize( ); u++ )
                         len -= L[ pairs[i].Left(u) ];
                    for ( int u = 0; u < pairs[i].RightSize( ); u++ )
                         len -= L[ pairs[i].Right(u) ];
                    double d = len / pairs[i].Dev( );
                    static vec<Bool> covered;
                    covered.resize_and_set( c.size( ), False );
                    for ( int k = 0; k < pairsx.isize( ); k++ )
                    {    const pp_pair& q = pairsx[k];
                         static vec<int> loffsets, roffsets;
                         GetOverlaps( c, q.Left( ), loffsets );
                         GetOverlaps( c, q.Right( ), roffsets );
                         for ( int ul = 0; ul < loffsets.isize( ); ul++ )
                         {    int lo = loffsets[ul];
                              for ( int ur = 0; ur < roffsets.isize( ); ur++ )
                              {    int ro = roffsets[ur];
                                   int g = ImpliedGap( c, q, L, lo, ro );
                                   if ( Abs( g - q.Gap( ) ) > dmult * q.Dev( ) )
                                        continue;
                                   int start = Max( 0, lo );
                                   int stop = Min( c.isize( ), lo + q.LeftSize( ) );
                                   for ( int r = start; r < stop; r++ )
                                        covered[r] = True;
                                   start = Max( 0, ro );
                                   stop = Min( c.isize( ), ro + q.RightSize( ) );
                                   for ( int r = start; r < stop; r++ )
                                        covered[r] = True;    }    }    }
                    if ( Sum(covered) == c.isize( ) ) 
                    {    ++valids;    
                         closures[i].push_back(c);    
                         devs[i].push_back(d);    }    }
               */

               if ( verbosity >= 1 )
                    cout << valids << " closures are valid" << endl;    }    }    

     if ( verbosity >= 1 )
          cout << "\n" << TimeSince(eclock) << " used" << endl;    }

void ClosePairedPairs( const vec<pp_pair>& pairs, const vec<Bool>& pairs_to_close,
     const vec<int>& L, const double dmult, vec< vec<pp_closure> >& closures, 
     vec< vec<double> >& devs, int max_ext, vec<Bool>& fail, int max_processed, 
     int max_unprocessed, int verbosity, int simple_walk_verbosity, int max_opens,
     int max_nodes, Bool create_closures_if_fail, const int max_closures )
{    fail.resize_and_set( pairs.size( ), False );
     closures.clear( ), devs.clear( );
     closures.resize( pairs.size( ) );
     devs.resize( pairs.size( ) );

     static vec<pp_read> reads;
     reads.clear( );
     for ( int i = 0; i < pairs.isize( ); i++ )
          reads.push_back( pairs[i].Left( ), pairs[i].Right( ) );
     UniqueSort(reads);
     vec< vec<int> > readsp;
     for ( int u = 0; u < pairs.isize( ); u++ )
          readsp.push_back( pairs[u].Left( ), pairs[u].Right( ) );
     vec_overlap<int> over(readsp);

     for ( int i = 0; i < pairs.isize( ); i++ )
     {    if ( !pairs_to_close[i] ) continue;
          if ( verbosity >= 1 )
          {    cout << "\nclosing pair " << i << " of " << pairs.size( )
                    << " = " << pairs[i] << endl;    }
          pp_mpair p( pairs[i] );
          Bool wfail;
          SimpleWalkRight( pairs[i], reads, L, dmult, closures[i], max_opens, 
               create_closures_if_fail, wfail, simple_walk_verbosity,
               False, max_nodes );
          if ( max_closures > 0 && closures[i].isize( ) > max_closures )
          {    if ( verbosity >= 1 ) cout << "Too many closures to keep.\n";
               fail[i] = True;
               closures[i].resize(0);
               continue;    }
          if (wfail) fail[i] = True;

          // Validate closures.  Note that we're not checking
          // for stretching, but should be.

          double vclock = WallClockTime( );
          ValidateClosures( p, closures[i], pairs, over, readsp, L, dmult, 
               verbosity );
          if ( verbosity >= 1 )
               cout << TimeSince(vclock) << " used validating closures" << endl;

          // Report results.

          report:
          if ( verbosity == 1 )
               cout << "found " << closures[i].size( ) << " validated closures\n";
          if ( verbosity >= 2 )
               cout << "found " << closures[i].size( ) << " validated closures:\n";
          for ( int j = 0; j < closures[i].isize( ); j++ )
          {    double len = -pairs[i].Gap( );
               for ( int u = 0; u < closures[i][j].isize( ); u++ )
                    len += L[ closures[i][j][u] ];
               for ( int u = 0; u < pairs[i].LeftSize( ); u++ )
                    len -= L[ pairs[i].Left(u) ];
               for ( int u = 0; u < pairs[i].RightSize( ); u++ )
                    len -= L[ pairs[i].Right(u) ];
               double d = len / pairs[i].Dev( );
               devs[i].push_back(d);
               if ( verbosity >= 2 )
               {    cout << closures[i][j] << " (off by ";
                    RightPrecisionOut( cout, d, 1 );
                    cout << " devs)" << endl;    }    }    }    }

Bool operator==( const pp_mpair& p1, const pp_mpair& p2 )
{    return p1.Left( ) == p2.Left( ) && p1.Right( ) == p2.Right( )
          && p1.Gap( ) == p2.Gap( ) && p1.Dev( ) == p2.Dev( )
          && p1.LeftMark( ) == p2.LeftMark( )
          && p1.RightMark( ) == p2.RightMark( );    }

Bool operator<( const pp_mpair& p1, const pp_mpair& p2 )
{    if ( p1.Left( ) < p2.Left( ) ) return True;
     if ( p1.Left( ) > p2.Left( ) ) return False;
     if ( p1.Right( ) < p2.Right( ) ) return True;
     if ( p1.Right( ) > p2.Right( ) ) return False;
     if ( p1.Gap( ) < p2.Gap( ) ) return True;
     if ( p1.Gap( ) > p2.Gap( ) ) return False;
     if ( p1.Dev( ) < p2.Dev( ) ) return True;
     if ( p1.Dev( ) > p2.Dev( ) ) return False;
     if ( p1.LeftMark( ) < p2.LeftMark( ) ) return True;
     if ( p1.LeftMark( ) > p2.LeftMark( ) ) return False;
     if ( p1.RightMark( ) < p2.RightMark( ) ) return True;
     if ( p1.RightMark( ) > p2.RightMark( ) ) return False;
     return False;    }

ostream& operator<<( ostream& o, const istring& x )
{    int count = 0;
     for ( int i = 0; i < x.isize( ); i++ )
     {    if ( i > 0 ) 
          {    o << ".";
               ++count;    }
          if ( count >= 80 && istring::s_displayMode == istring::SPLIT )
          {    o << "\n";
               count = 0;    }
          o << BaseAlpha( x[i] );
          count += BaseAlpha( x[i] ).size( );    }
     return o;    }

void BinaryWrite( int fd, const pp_pair& p )
{    BinaryWrite( fd, p.Left( ) );
     BinaryWrite( fd, p.Right( ) );
     WriteBytes( fd, &p.gap_, sizeof(p.gap_) );
     WriteBytes( fd, &p.dev_, sizeof(p.dev_) );    }

void BinaryRead( int fd, pp_pair& p )
{    BinaryRead( fd, p.LeftMutable( ) );
     BinaryRead( fd, p.RightMutable( ) );
     ReadBytes( fd, &p.gap_, sizeof(p.gap_) );
     ReadBytes( fd, &p.dev_, sizeof(p.dev_) );    }

ostream& operator<<( ostream& o, const pp_pair& p )
{    o << p.Left( ) << " --- ";
     if ( p.Dev( ) == 0 && double( int( floor( p.Gap( ) ) ) ) == p.Gap( ) )
          o << int( floor( p.Gap( ) ) );
     else
     {    RightPrecisionOut( o, p.Gap( ), 2 );
          o << " +/- ";
          RightPrecisionOut( o, p.Dev( ), 2 );    }
     return o << " --> " << p.Right( );    }

void pp_pair::Print( ostream& o, const vec<int>& L ) const
{    if ( !IsClosed( *this, L ) ) o << *this;
     else o << this->Left( ) << " (" << -int( floor( this->Gap( ) ) ) << ")";    }

ostream& operator<<( ostream& o, const pp_mpair& p )
{    for ( int i = 0; i < p.Left( ).isize( ); i++ )
     {    if ( i == p.LeftMark( ) ) o << "|";
          else if ( i > 0 ) o << ".";
          o << BaseAlpha( p.Left(i) );    }
     if ( p.LeftMark( ) == p.LeftSize( ) ) o << "|";
     o << " --- ";
     RightPrecisionOut( o, p.Gap( ), 2 );
     o << " +/- ";
     RightPrecisionOut( o, p.Dev( ), 2 );
     o << " --> ";
     for ( int i = 0; i < p.Right( ).isize( ); i++ )
     {    if ( i == p.RightMark( ) ) o << "|";
          else if ( i > 0 ) o << ".";
          o << BaseAlpha( p.Right(i) );    }
     if ( p.RightMark( ) == p.RightSize( ) ) o << "|";
     return o;    }

void ExtendPairedPairs( const vec<pp_pair>& pairs, const vec<int>& L,     
     const double dmult, vec<pp_pair>& extensions )
{    extensions.resize( pairs.size( ) );
     for ( int i = 0; i < pairs.isize( ); i++ )
     {    Bool infinite_loop = False;
          pp_mpair p( pairs[i] );
          cout << "\n";
          PRINT3( i, pairs.size( ), p );
          Bool try_close = True;
          restart:
          while(1)
          {    if ( try_close && ReadyToClose( p, dmult ) )
               {    static vec<pp_closure> cl;
                    Bool trim = False;
                    GetClosures( p, L, dmult, cl, False );
                    if ( cl.size( ) == 1 )
                    {    double g = 0.0;
                         for ( int j = 0; j < cl[0].isize( ); j++ )
                              g -= L[ cl[0][j] ];
                         if ( cl[0].isize( ) > p.LeftSize( )
                              || cl[0].isize( ) > p.RightSize( )
                              || p.Gap( ) != g || p.Dev( ) != 0.0 )
                         {    p = pp_mpair( pp_pair( cl[0], cl[0], g, 0.0 ) );
                              continue;    }    }    }
               static vec<pp_mpair> pnew;
               FindExtensions( p, pairs, L, dmult, pnew );
               static vec<int> leftleft, leftright, rightleft, rightright; 
               leftleft.clear( ), leftright.clear( );
               rightleft.clear( ), rightright.clear( );
               for ( int j = 0; j < pnew.isize( ); j++ )
               {    const pp_mpair& y = pnew[j];
                    if ( y.LeftMark( ) > 0 )
                    {    leftleft.push_back( y.Left( y.LeftMark( ) - 1 ) );    }
                    if ( y.LeftSize( ) - y.LeftMark( ) > p.LeftSize( ) )
                    {    leftright.push_back(
                              y.Left( p.LeftSize( ) + y.LeftMark( ) ) );    }
                    if ( y.RightMark( ) < y.RightSize( ) )
                    {    rightright.push_back( y.Right( y.RightMark( ) ) );    }
                    if ( y.RightMark( ) > p.RightSize( ) )
                    {    rightleft.push_back( y.Right( 
                              y.RightMark( ) - p.RightSize( ) - 1 ) );    }    }
               UniqueSort(leftleft), UniqueSort(leftright);
               UniqueSort(rightleft), UniqueSort(rightright);
               if ( leftleft.size( ) != 1 && leftright.size( ) != 1 
                    && rightleft.size( ) != 1 && rightright.size( ) != 1 )
               {    break;    }
               pp_read left = p.Left( ), right = p.Right( );
               double g = p.Gap( );
               if ( leftleft.size( ) == 1 ) left.push_front( leftleft[0] );
               if ( leftright.size( ) == 1 ) 
               {    left.push_back( leftright[0] );
                    g -= L[ leftright[0] ];    }
               if ( rightleft.size( ) == 1 ) 
               {    right.push_front( rightleft[0] );
                    g -= L[ rightleft[0] ];    }
               if ( rightright.size( ) == 1 ) right.push_back( rightright[0] );
               p = pp_mpair( pp_pair( left, right, g, p.Dev( ) ) );
               if ( p.LeftSize( ) + p.RightSize( ) > 100 )
               {    p = pairs[i];
                    if (try_close)
                    {    try_close = False;
                         goto restart;    }
                    infinite_loop = True;
                    cout << "maybe in infinite loop\n";
                    break;    }    }
          extensions[i] = p;
          PRINT(p);    }    }

// ClosePairedPairsDirectUnique: if the easy graph reveals a unique closure for a
// given pair, close it.  This does not do any validation, but perhaps it should.
// The complexity of the allowed computation is controlled by max_partials.
// Perhaps this should be combined with BridgeGaps.
//
// Also, if a pair has a small number of closures, all of the same length, adjust
// its gap and deviation accordingly.
//
// Return number of changes.

int ClosePairedPairsDirectUnique( vec<pp_pair>& ppp, const HyperKmerPath& h,
     int max_loops, const vec<Bool>& remove )
{    vec<int> L;
     vec<int> to_left, to_right;
     for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
          L.push_back( h.EdgeObject(i).KmerCount( ) );
     digraphE<int> G( h, L );
     double dmult = 4.0;
     G.ToLeft(to_left), G.ToRight(to_right);
     int changes = 0;
     const int max_closures = 5;
     for ( int i = 0; i < ppp.isize( ); i++ )
     {    if ( remove[i] || ppp[i].Dev( ) == 0 ) continue;
          static vec< vec<int> > paths;
          int v = to_right[ ppp[i].Left( ).back( ) ]; 
          int w = to_left[ ppp[i].Right( ).front( ) ];
          int L1 = int( ceil( ppp[i].Gap( ) - dmult * ppp[i].Dev( ) ) );
          int L2 = int( floor( ppp[i].Gap( ) + dmult * ppp[i].Dev( ) ) );
          static vec<pp_closure> closures;
          closures.clear( );
          if ( L1 < 0 )
          {    GetClosures( ppp[i], L, dmult, closures, False );
               L1 = 0;    }
          if ( G.AllPathsLengthRange( v, w, L1, L2, to_right, paths, max_closures, 
               max_loops ) )
          {    for ( int j = 0; j < paths.isize( ); j++ )
               {    pp_closure c = ppp[i].Left( );
                    c.append( paths[j] );
                    c.append( ppp[i].Right( ) );
                    closures.push_back(c);    }
               if ( closures.empty( ) ) continue;
               pp_read c = closures[0];
               int len = 0;
               for ( int j = 0; j < c.isize( ); j++ )
                    len += L[ c[j] ];
               if ( closures.solo( ) )
               {    ppp[i] = pp_pair( c, c, -double(len), 0.0 );    
                    ++changes;    }
               else
               {    int x;
                    for ( x = 0; x < closures.isize( ); x++ )
                    {    int lenx = 0;
                         for ( int j = 0; j < closures[x].isize( ); j++ )
                              lenx += L[ closures[x][j] ];
                         if ( lenx != len ) break;    }
                    if ( x == closures.isize( ) )
                    {    int lr = 0;
                         for ( int j = 0; j < ppp[i].LeftSize( ); j++ )
                              lr += L[ ppp[i].Left(j) ];
                         for ( int j = 0; j < ppp[i].RightSize( ); j++ )
                              lr += L[ ppp[i].Right(j) ];
                         ppp[i] = pp_pair( ppp[i].Left( ), ppp[i].Right( ),
                              double(len-lr), 0.0 );    
                         ++changes;    }    }    }    }
     return changes;    }

Bool IsClosed( const pp_pair& p, const vec<int>& L )
{    if ( p.Dev( ) != 0 ) return False;
     if ( p.Left( ) != p.Right( ) ) return False;
     int len = 0;
     for ( int i = 0; i < p.LeftSize( ); i++ )
          len += L[ p.Left(i) ];
     return p.Gap( ) == -double(len);    }
