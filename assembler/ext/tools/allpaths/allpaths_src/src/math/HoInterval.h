///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef HO_INTERVAL_H
#define HO_INTERVAL_H

#include "Vec.h"
#include "math/Functions.h"
#include "feudal/MasterVec.h"
#include "feudal/SerfVec.h"

/// Class: ho_interval
/// A half-open interval [a, b).
class ho_interval {

     public:

     ho_interval( ) { }
     ho_interval( int start, int stop ) : start_(start), stop_(stop) 
     {    ForceAssertGe( static_cast<long>(stop)-start, 0 );    }

     int Start( ) const { return start_; }
     int Stop( ) const { return stop_; }

     int Length( ) const 
     {    ForceAssertLe( start_, stop_ );
          return stop_ - start_;    }

     void Set( int start, int stop )
     {    start_ = start;
          stop_ = stop;    }
     void SetStart( int start )
     {    start_ = start;    }
     void SetStop( int stop )
     {    stop_ = stop;    }
     void AddToStart( int delta )
     {    start_ += delta;    }
     void AddToStop( int delta )
     {    stop_ += delta;    }
     void Shift( int delta )
     {    start_ += delta;
          stop_ += delta;    }
     // Lengthen interval by lengthen on each side (if negative,
     // shorten).  Stay in [0,maxStop), and if length becomes <=0 make
     // start = old stop.
     void Lengthen( int lengthen, int maxStop )
     {    if (0!=lengthen)
          {    start_ = min(max(0, start_ - lengthen), stop_);
               stop_ = max(min(stop_ + lengthen, maxStop), start_);    }    }
     bool Contains(int i) const { return i >= start_ && i < stop_; }

     /// If the two overlap, merge into this and return true.
     bool Merge(const ho_interval & o);

     friend Bool operator==( const ho_interval& i1, const ho_interval& i2 )
     {    return i1.start_ == i2.start_ && i1.stop_ == i2.stop_;    }

     friend Bool operator!=( const ho_interval& i1, const ho_interval& i2 )
     {    return !( i1 == i2 );    }

     friend Bool operator<( const ho_interval& i1, const ho_interval& i2 )
     {    return i1.Start( ) < i2.Start( ) 
               || ( i1.Start( ) == i2.Start( ) && i1.Stop( ) < i2.Stop( ) );    }

     friend Bool operator>( const ho_interval& i2, const ho_interval& i1 )
     {    return i1.Start( ) < i2.Start( ) 
               || ( i1.Start( ) == i2.Start( ) && i1.Stop( ) < i2.Stop( ) );    }

     friend ho_interval operator-( const ho_interval& h, int x )
     {    return ho_interval( h.Start( ) - x, h.Stop( ) - x );    }

     friend ho_interval operator+( const ho_interval& h, int x )
     {    return ho_interval( h.Start( ) + x, h.Stop( ) + x );    }

     friend ostream&operator<<( ostream& out, const ho_interval& i )
     {    return out << i.Start( ) << "-" << i.Stop( );    }

     protected:

     // TODO: potentially dangerous truncation of indexes
     int start_, stop_;

};

TRIVIALLY_SERIALIZABLE(ho_interval);
typedef SerfVec<ho_interval> HOIntervalVec;
typedef MasterVec<HOIntervalVec> VecHOIntervalVec;

///Compare two ho_intervals by length.
struct SmallerLength: public binary_function<ho_interval, ho_interval, bool> {
  bool operator()(const ho_interval & l, const ho_interval & r) {
    return (l.Length() < r.Length());
  }
};

inline longlong Sum( const vec<ho_interval>& v )
{    longlong sum = 0;
     for ( size_t i = 0; i < v.size(); i++ )
          sum += v[i].Length( );
     return sum;    }

inline Bool Member( const ho_interval& x, int k )
{    return k >= x.Start( ) && k < x.Stop( );    }

inline Bool Member( const vec<ho_interval>& v, int k )
{    for ( size_t i = 0; i < v.size(); i++ )
          if ( Member( v[i], k ) ) return True;
     return False;    }

/// Find out whether int k is included in a sorted vector.
/// Faster than Member: log instead of linear.
/// maxSize is the Length() of the largest interval in the vector. It should be
/// calculated ahead of time, because if we calculate it inline the 
/// function becomes linear again!

inline int PositionSorted( const vec<ho_interval>& v, int k, int maxSize ) {
  ho_interval toFind(k+1,k+2);
  int end = upper_bound(v.begin(), v.end(), toFind) - v.begin();
  end = min(end, (int)v.size()-1);
  int start = k - maxSize;
  for (int i=end; i >= 0; --i) {
    if (v[i].Contains(k)) return i;
    if (v[i].Start() < start) break;
  }
  return -1;
}

inline bool MemberSorted( const vec<ho_interval>& v, int k, int maxSize )
{    return PositionSorted( v, k, maxSize ) >= 0;    }

///Return positive distance if they don't overlap, 0 otherwise.
inline int Distance( const ho_interval& x, const ho_interval& y )
{    if ( x.Stop( ) < y.Start( ) ) return y.Start( ) - x.Stop( );
     if ( y.Stop( ) < x.Start( ) ) return x.Start( ) - y.Stop( );
     return 0;    }

/// Function: Distance
/// Return the distance between two half-open intervals (from end of one to start
/// of other) if they don't overlap, 0 otherwise.
/// This method pretends that [x1,x2) and [y1,y2) are HoIntervals.
inline int Distance( int x1, int x2, int y1, int y2 )
{    if ( x2 < y1 ) return y1 - x2;
     if ( y2 < x1 ) return x1 - y2;
     return 0;    } 

/// Return length of overlap, 0 if there is none.
int Overlap( const ho_interval& x, const ho_interval& y );

/// Meets( x, y ): return true if two ho_intervals x, y overlap.

inline Bool Meets( const ho_interval& x, const ho_interval& y )
{    return Overlap( x, y ) > 0;    }

/// Overlap.  Suppose that i2 is an ordered list of disjoint intervals.
/// Compute the sum of Overlap( i1, i2[x] ) as x varies.

int Overlap( const ho_interval& i1, const vec<ho_interval>& i2 );

/// OverlapIndices( x, v, L, I ): for a ho_interval x and a sorted 
/// vec<ho_interval> v, for which the maximum length is L, find all y in v such that 
/// x meets y, and return the corresponding indices as I.

void OverlapIndices( const ho_interval& x, const vec<ho_interval>& v, const int L,
     vec<int>& I );

/// operator^: take intersection of two ho_intervals.  Note that the answer
/// isn't proper if the intervals don't overlap.

inline ho_interval operator^( const ho_interval& i1, const ho_interval& i2 )
{    return ho_interval( 
          Max( i1.Start( ), i2.Start( ) ),
          Min( i1.Stop( ), i2.Stop( ) ) );    }

inline ho_interval Span( const ho_interval& i1, const ho_interval& i2 )
{    return ho_interval( 
          Min( i1.Start( ), i2.Start( ) ),
          Max( i1.Stop( ), i2.Stop( ) ) );    }

/// Subset: determine if one ho_interval is contained in another.

inline Bool Subset( const ho_interval& i1, const ho_interval& i2 )
{    return i1.Start( ) >= i2.Start( ) && i1.Stop( ) <= i2.Stop( );    }

/// Subset: determine if one ho_interval is contained in the union of some others.

Bool Subset( const ho_interval& i1, const vec<ho_interval>& i2 );

/// CondenseIntervals: given a set of half-open intervals from [0,n), convert it
/// to a sorted list of half-open intervals, whose disjoint union is [0,n), and to
/// each such interval pair the coverage by the original intervals.

void CondenseIntervals( int n, const vec<ho_interval>& in, 
     vec< pair<ho_interval, int> >& out );

/// TotalCovered: return the total coverage of a list of half-open intervals of
/// of nonnegative integers.

int TotalCovered( const vec<ho_interval>& h );

/// ExtractGivenCoverage: given a list of half-open intervals from [0,n), return an
/// ordered list of disjoint intervals on which the coverage by the given list is
/// at least c.  The case where in = out is allowed.

void ExtractGivenCoverage( int n, int c, const vec<ho_interval>& in,
     vec<ho_interval>& out );

/// Uncovered: given a list of half-open intervals from [0,n), return ordered
/// list of disjoint intervals whose union is the uncovered part of [0,n).

void Uncovered( int n, const vec<ho_interval>& in, vec<ho_interval>& un );

/// RemoveNearDuplicates: given a list of intervals, remove those which nearly
/// overlap each other (leaving only a single copy).  Two intervals nearly 
/// overlap each other if their left endpoints differ by at most max_diff, 
/// and so do their right endpoints.  Exactly which intervals are removed 
/// is not well-defined.

void RemoveNearDuplicates( vec<ho_interval>& vi, int max_diff );

int Span( const vec<ho_interval>& v );










//============== class HoIntervalWithId =============================//



///This class keeps track of where the interval came from with an int id.
///Useful when sorting intervals and wanting to refer back to their origin.
class HoIntervalWithId: public ho_interval {
public:
  int id; // TODO: potentially dangerous truncation of index
  HoIntervalWithId( int start=0, int stop=1, int id=0 ) : 
    ho_interval(start, stop), id(id) {}

  ///Read it in from three-column format.
  friend istream & operator>>(istream & is, HoIntervalWithId & i) {
    is >> i.id >> i.start_ >> i.stop_;
    return is;
  }

  bool Contains(const HoIntervalWithId & o) const {
    return id == o.id && start_ <= o.start_ && stop_ >= o.stop_;
  }

  bool Merge(const HoIntervalWithId & o) {
    if (id != o.id) return false;
    else return ho_interval::Merge(o);
  }

  friend Bool operator==(const HoIntervalWithId & a, const HoIntervalWithId & b) {
    return ( a.id == b.id ) && ( a.start_ == b.start_ ) && ( a.stop_ == b.stop_ );
  }

  HoIntervalWithId & operator += ( int x ) {
    start_ += x ;
    stop_  += x ;
    return *this;
  }

  HoIntervalWithId & operator -= ( int x ) {
    start_ -= x ;
    stop_  -= x ;
    return *this;
  }

  friend HoIntervalWithId operator-( const HoIntervalWithId & h, int x )
  {    return HoIntervalWithId( h.Start( ) - x, h.Stop( ) - x, h.id );    }

  friend HoIntervalWithId operator+( const HoIntervalWithId& h, int x )
  {    return HoIntervalWithId( h.Start( ) + x, h.Stop( ) + x, h.id );    }


};

/// operator^: take intersection of two ho_intervals.  Will return an empty interval (start==stop) if the operands
/// don't overlap; in this case the actual values of start and stop are not defined and should not be relied upon

inline HoIntervalWithId operator^( const HoIntervalWithId & i1, const HoIntervalWithId & i2 ) {  
    if ( i1.id != i2.id ) return HoIntervalWithId(0,0,i1.id);
    int l = Max( i1.Start( ), i2.Start( ) );
    int r = Min( i1.Stop( ), i2.Stop( ) ) ;    
    if ( l < r ) return HoIntervalWithId(l,r,i1.id);
    return HoIntervalWithId(0,0,i1.id);
}


///Order by id and then by start.
inline 
bool LessById(const HoIntervalWithId & l, const HoIntervalWithId & r) { 
  if (l.id == r.id) {
    return l < r; 
  } 
  return l.id < r.id;
}

///Order by id only.
inline 
bool LessByIdOnly(const HoIntervalWithId & l, const HoIntervalWithId & r) { 
  return l.id < r.id;
}

///Read in a HoIntervalWithId as three columns.
inline 
ostream & operator<< (ostream & os, const HoIntervalWithId & i) {
  return os << i.id << "\t" << i.Start() << "\t" << i.Stop();
}

/// Take in a sorted vector, merge together all overlapping intervals.
void Merge(vec<HoIntervalWithId> & v);

/// Return length of overlap, 0 if there is none.
int Overlap( const HoIntervalWithId & x, const HoIntervalWithId & y );

/// Meets( x, y ): return true if two intervals x, y overlap.
inline Bool Meets( const HoIntervalWithId & x, const HoIntervalWithId & y ) {    
  return Overlap( x, y ) > 0;    
}

/// Meets( x, y ): return true if any two intervals from within x, y overlap.
Bool Meets( const vec<HoIntervalWithId> & x, const vec<HoIntervalWithId> & y );

///  Returns <true> if explicitly specified position <contig>:<pos> is to the left from the interval
/// start (in generalized sense: if it is at lower contig, it is less too!)
inline bool Less( unsigned int contig, unsigned int pos, const HoIntervalWithId & h) {
  return   ( contig < (unsigned int)h.id || ( contig == (unsigned int)h.id && pos < (unsigned int)h.Start() ) );
}

///  Returns <true> if explicitly specified position <contig>:<pos> is to the right from the interval
/// end (in generalized sense: if it is at higher contig, it is greater too!)
inline bool Greater( unsigned int contig, unsigned int pos, const HoIntervalWithId & h) {
  return   ( contig > (unsigned int)h.id || ( contig == (unsigned int)h.id && pos >= (unsigned int)h.Stop() ) );
}

///  Returns <true> if explicitly specified position <contig>:<pos> is within the interval's boundaries
inline bool Inside( unsigned int contig, unsigned int pos, const HoIntervalWithId & h) {
  return   ( contig == (unsigned int)h.id && pos >= (unsigned int)h.Start() && pos < (unsigned int)h.Stop() );
}

namespace std {
template<>
struct less< vec<HoIntervalWithId> > : 
  public binary_function< const vec<HoIntervalWithId> &,  const vec<HoIntervalWithId> &, bool> {
  bool operator() (const vec<HoIntervalWithId> & a,  const vec<HoIntervalWithId> & b) const {

    unsigned int min_size = min(a.size(),b.size());

    for ( unsigned int i = 0 ; i < min_size ; i++ ) {
      if (LessById(a[i],b[i])) return true;
      if (LessById(b[i],a[i])) return false;
    }
    // we have compared intervals up to min_size by now;
    // if we got here, all those intervals were the same

    // so now we check the sizes; the shorter vector is "less"
    if ( a.size() < b.size() ) return true;
    if ( a.size() > b.size() ) return false;

    // all the intervals up to min_size were the same, and the sizes
    // of a and b are equal, so that the vectors are exactly equal:
    return false;
  }
};
}


/*
inline bool MemberSorted(const vec<HoIntervalWithId> & intervals, int id,
			 int pos) {
  typedef vec<HoIntervalWithId>::iterator viter;
  HoIntervalWithId contig(0, 0, id);
  pair<viter,viter> contigRange = 
    equal_range(intervals.begin(), intervals.end(), contig, LessByIdOnly);
*/
  
  


#endif
