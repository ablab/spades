///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "math/Functions.h"
#include "math/HoInterval.h"
#include "STLExtensions.h"
#include "feudal/OuterVecDefs.h"
#include "feudal/SmallVecDefs.h"

template class SmallVec< ho_interval, MempoolAllocator<ho_interval> >;
template class OuterVec<HOIntervalVec>;

bool ho_interval::Merge(const ho_interval & o) {
  if (!Meets(*this,o)) return false;
  *this = Span(*this, o);
  return true;
}

int Overlap( const ho_interval& x, const ho_interval& y ) {    
  return Max( 0, 
          Min( x.Stop( ), y.Stop( ) ) - Max( x.Start( ), y.Start( ) ) );    
}


/// Return length of overlap, 0 if there is none.
int Overlap( const HoIntervalWithId & x, const HoIntervalWithId & y ) {
  if ( x.id != y.id ) return 0;
  return Max( 0, 
          Min( x.Stop( ), y.Stop( ) ) - Max( x.Start( ), y.Start( ) ) );    
}

/// Meets( x, y ): return true if any two intervals from within x, y overlap.
Bool Meets( const vec<HoIntervalWithId> & x, const vec<HoIntervalWithId> & y ) {
  for ( unsigned int i = 0 ; i < x.size() ; i++ ) {
    for ( unsigned int j = 0 ; j < y.size() ; j++ ) {
      if ( Meets(x[i],y[j]) ) return true;
    }
  }
  return false;
}


// CondenseIntervals: given a set of half-open intervals from [0,n), convert it
// to a sorted list of half-open intervals, whose disjoint union is [0,n), and to
// each such interval pair the coverage by the original intervals.

void CondenseIntervals( int n, const vec<ho_interval>& in, 
     vec< pair<ho_interval, int> >& out )
{    vec< pair<int,int> > start_stop;
     for ( int i = 0; i < (int) in.size( ); i++ )
     {    start_stop.push_back( make_pair( in[i].Start( ), +1 ) );
          start_stop.push_back( make_pair( in[i].Stop( ), -1 ) );    }
     sort( start_stop.begin( ), start_stop.end( ) );
     out.clear( );
     if ( start_stop.size( ) == 0 )
     {    out.push_back( make_pair( ho_interval( 0, n ), 0 ) );
          return;    }
     if ( start_stop[0].first != 0 )
          out.push_back( make_pair( ho_interval( 0, start_stop[0].first ), 0 ) );
     int coverage = 0;
     for ( int i = 0; i < (int) start_stop.size( ) - 1; i++ )
     {    coverage += start_stop[i].second;
          if ( start_stop[i].first != start_stop[i+1].first )
               out.push_back( make_pair( 
                    ho_interval( start_stop[i].first, start_stop[i+1].first ), 
                         coverage ) );    }
     coverage += start_stop.back( ).second;
     if ( start_stop.back( ).first < n )
          out.push_back( make_pair( ho_interval( start_stop.back( ).first, n ),
               coverage ) );    }

void ExtractGivenCoverage( int n, int c, const vec<ho_interval>& in,
     vec<ho_interval>& out )
{    vec< pair<ho_interval, int> > condensed;
     CondenseIntervals( n, in, condensed );
     out.clear( );
     for ( int i = 0; i < (int) condensed.size( ); i++ )
          if ( condensed[i].second >= c )
          {    int j;
               for ( j = i + 1; j < (int) condensed.size( ); j++ )
                    if ( condensed[j].second < c ) break;
               out.push_back( ho_interval( condensed[i].first.Start( ), 
                    condensed[j-1].first.Stop( ) ) );    
               i = j;    }    }

// RemoveNearDuplicates: given a list of intervals, remove those which nearly
// overlap each other (leaving only a single copy).  Two intervals nearly overlap
// each other if their left endpoints differ by at most max_diff, and so do their
// right endpoints.  Exactly which intervals are removed is not well-defined.

void RemoveNearDuplicates( vec<ho_interval>& vi, int max_diff )
{    sort( vi.begin( ), vi.end( ) );
     vec<Bool> keep( vi.size( ), True );
     for ( int i = 0; i < (int) vi.size( ); i++ )
     {    if ( !keep[i] ) continue;
          int j;
          for ( j = i + 1; j < (int) vi.size( ); j++ )
               if ( vi[j].Start( ) - vi[i].Start( ) > max_diff ) break;
          for ( int k = i + 1; k < j; k++ )
               if ( Abs( vi[k].Stop( ) - vi[i].Stop( ) ) <= max_diff )
                    keep[k] = False;    }
     int count = 0;
     for ( int i = 0; i < (int) vi.size( ); i++ )
     {    if ( !keep[i] ) continue;
          if ( i != count ) vi[count] = vi[i];
          ++count;    }
     vi.resize(count);    }

Bool Subset( const ho_interval& i1, const vec<ho_interval>& i2 )
{    if ( i2.empty( ) ) return False;
     int m = i2[0].Start( ), M = i2[0].Stop( );
     for ( size_t k = 1; k < i2.size(); k++ )
     {    m = Min( m, i2[k].Start( ) );
          M = Max( M, i2[k].Stop( ) );    }
     vec<ho_interval> i2copy = i2;
     for ( size_t k = 0; k < i2.size(); k++ )
          i2copy[k].Shift(-m);
     ho_interval i1copy;
     i1copy = i1;
     i1copy.Shift(-m);
     vec<ho_interval> out;
     ExtractGivenCoverage( M - m, 1, i2copy, out );
     for ( size_t k = 0; k < out.size(); k++ )
          if ( Subset( i1copy, out[k] ) ) return True;
     return False;    }

// Overlap is implemented inefficiently.  It should do a binary search.

int Overlap( const ho_interval& i1, const vec<ho_interval>& i2 )
{    int start = lower_bound( i2.begin( ), i2.end( ), i1 ) - i2.begin( ) - 1;
     int sum = 0;
     size_t k = start > 0 ? start : 0;
     for ( ; k < i2.size(); k++ )
     {    if ( i2[k].Start( ) >= i1.Stop( ) ) break;
          sum += Overlap( i1, i2[k] );    }
     return sum;    }

int TotalCovered( const vec<ho_interval>& h )
{    vec< pair<ho_interval, int> > out;
     int M = 0;
     for ( int i = 0; i < (int) h.size( ); i++ )
          M = Max( M, h[i].Stop( ) );
     CondenseIntervals( M, h, out );
     int answer = 0;
     for ( int i = 0; i < (int) out.size( ); i++ )
          if ( out[i].second > 0 ) answer += out[i].first.Length( );
     return answer;    }

int Span( const vec<ho_interval>& v )
{    ForceAssert( v.nonempty( ) );
     vec<int> left, right;
     for ( size_t i = 0; i < v.size(); i++ )
     {    left.push_back( v[i].Start( ) );
          right.push_back( v[i].Stop( ) );    }
     return Max(right) - Min(left);    }

void Uncovered( int n, const vec<ho_interval>& in, vec<ho_interval>& un )
{    vec<ho_interval> cov;
     ExtractGivenCoverage( n, 1, in, cov );
     un.clear( );
     if ( cov.empty( ) ) 
     {    un.push_back( ho_interval( 0, n ) );
          return;    }
     if ( cov[0].Start( ) > 0 ) un.push_back( ho_interval( 0, cov[0].Start( ) ) );
     for ( size_t j = 0; j+1 < cov.size(); j++ )
          un.push_back( ho_interval( cov[j].Stop( ), cov[j+1].Start( ) ) );
     if ( cov.back( ).Stop( ) < n )
          un.push_back( ho_interval( cov.back( ).Stop( ), n ) );    }

/// OverlapIndices( x, v, L, I ): for a ho_interval x and a sorted 
/// vec<ho_interval> v, for which the maximum length is L, find all y in v such that
/// x meets y, and return the corresponding indices as I.

void OverlapIndices( const ho_interval& x, const vec<ho_interval>& v, const int L,
     vec<int>& I )
{    I.clear( );
  ho_interval h( Min( x.Stop( ) - 1, x.Stop( ) ), x.Stop( ) );
  int end = upper_bound( v.begin( ), v.end( ), h ) - v.begin( );
  int start = x.Start( ) - L;
  for ( int i = end - 1; i >= 0; --i ) 
  {    
    if ( Meets( x, v[i] ) ) I.push_back(i);
    if ( v[i].Start( ) <= start ) return;    
  }
}

void Merge(vec<HoIntervalWithId> & v) {
  Assert(is_sorted(v.begin(), v.end(), LessById));
  HoIntervalWithId BAD(-2,-1,-1);
  for (size_t i=1; i < v.size(); ++i) {
    if (v[i].Merge(v[i-1]) ) { v[i-1] = BAD; }
  }
  v.erase(remove(v.begin(), v.end(), BAD), v.end());
}
