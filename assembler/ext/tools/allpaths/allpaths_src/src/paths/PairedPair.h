/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2005) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// This file defines a function ClosePairedPairs that uses a kind of insert
// walking to find all closures of read pairs, provided that those read pairs
// have been turned into unipath sequences.  The function follows the 
// "paired pair" requirement, that to use a read pair to help close an insert,
// it must overlap on both ends.  This gains specificity and the expense of
// sensitivity.
// 
// What is given here is a "minimal" mathematical description of the problem.

/* Definition of the problem.  

1. First, a read (class pp_read) is a sequence of integers from [0,N), where 
N is specified.  Every element x of [0,N) is assigned an integer length L[x].
A read pair (pp_pair) is a pair of reads, together with a gap between them,
of size s +/- d.  A marked read pair (class pp_mpair) is a read pair, together 
with a position on its left read and a position on its right read.

2. We are given a collection S of read pairs.  The goal is to take each in turn,
and close its gap by filling in with sequence from other reads, according to
the following rules:
(a) the given read pair is first turned into a marked read pair X, where the 
marks are before the beginning of the left read, and after the end of the 
right read;
(b) iteratively, we then extend X by reaching into S and merging it with
a matching read pair P, yielding X' (which replaces X), e.g.
             1,7,4,3   -->  3,5,5       X  (gap = s +/- d)
                 4,3,6 -->      5,7,2   P  (gap = t +/- e)
             1,7,4,3,6 -->  3,5,5,7,2   X' (gap = (s - L[6]) +/- d)
where the marks on X are carried forward to X';
(c) X and P must match (align end to end) on both sides;
(d) such a merger is only allowed if the implied gap g1 on X' as computed from
X, e.g. s - L[6] is compatible with the gap g2 on X' as computed from
P, e.g. t - L[3] - L[5], in the sense that |g1-g2| <= dmult(d+e), where dmult 
is a fixed constant;
(e) if after a given merger, or initially, the gap g on X is negative by more
than E = dmult * d, a gap-closure is required: we find all possible ways of 
joining the left and right reads, requiring that the overlap (computed via L) 
is within E of g;
(f) after making such a closure, we trim off the parts of the sequences that
are to the left of the left mark point and to the right of the right mark
point;
(g) a closure (class pp_closure) is thus a sequence of integers from [0,N);
(h) subject to these rules, every possible closure is found.
*/

#ifndef PAIRED_PAIR_H
#define PAIRED_PAIR_H

#include "CoreTools.h"

class istring : public vec<int> {

     public:

     istring( ) { }
     explicit istring( const vec<int>& v ) : vec<int>(v) { }

     int Length( const vec<int>& L ) const;

     void Append( const istring& x )
     {    append(x);    }

 private:
     enum display_mode { ONE_LINE, SPLIT };
     static display_mode s_displayMode;

 public:

     // By default, istrings are outputted in multi-line format (splitting so that
     // no one line is too long).  You can override this with:
     //     istring::UseOneLineOutput( );
     static void UseOneLineOutput() { s_displayMode = ONE_LINE; }
     static void UseSplitLineOutput() { s_displayMode = SPLIT; }

     friend ostream& operator<<( ostream& o, const istring& x );

     friend void BinaryWrite( int fd, const istring& x )
     {    BinaryWrite( fd, (const vec<int>&) x );    }
     friend void BinaryRead( int fd, istring& x )
     {    BinaryRead( fd, (vec<int>&) x );    }

};

inline void BinaryWrite( int fd, const vec<istring>& v )
{    BinaryWriteComplex( fd, v );    }
inline void BinaryRead( int fd, vec<istring>& v )
{    BinaryReadComplex( fd, v );    }

typedef istring pp_read;

class pp_pair {

     public:

     pp_pair( ) { }
     pp_pair( const pp_read& left, const pp_read& right, double gap, double dev )
          : left_(left), right_(right), gap_(gap), dev_(dev) { }

     const pp_read& Left( ) const { return left_; }
     const pp_read& Right( ) const { return right_; }
     pp_read& LeftMutable( ) { return left_; }
     pp_read& RightMutable( ) { return right_; }
     int LeftSize( ) const { return left_.size( ); }
     int RightSize( ) const { return right_.size( ); }
     int Left( int i ) const { return left_[i]; }
     int Right( int i ) const { return right_[i]; }
     void SetLeft( const pp_read& p ) { left_ = p; }
     void SetRight( const pp_read& p ) { right_ = p; }
     
     double Gap( ) const { return gap_; }
     void SetGap( double g ) { gap_ = g; }
     double Dev( ) const { return dev_; }
     void SetDev( double d ) { dev_ = d; }

     void ReverseMe( )
     {    left_.ReverseMe( );
          right_.ReverseMe( );
          swap( left_, right_ );    }

     friend ostream& operator<<( ostream& o, const pp_pair& p );
     void Print( ostream& o, const vec<int>& L ) const;

     friend void BinaryWrite( int fd, const pp_pair& p );
     friend void BinaryRead( int fd, pp_pair& p );

     friend Bool operator<( const pp_pair& p1, const pp_pair& p2 )
     {    if ( p1.left_ < p2.left_ ) return True;
          if ( p1.left_ > p2.left_ ) return False;
          if ( p1.right_ < p2.right_ ) return True;
          if ( p1.right_ > p2.right_ ) return False;
          if ( p1.gap_ < p2.gap_ ) return True;
          if ( p1.gap_ > p2.gap_ ) return False;
          if ( p1.dev_ < p2.dev_ ) return True;
          return False;    }

     friend Bool operator==( const pp_pair& p1, const pp_pair& p2 )
     {    return p1.left_ == p2.left_ && p1.right_ == p2.right_
               && p1.gap_ == p2.gap_ && p1.dev_ == p2.dev_;    }
     friend Bool operator!=( const pp_pair& p1, const pp_pair& p2 )
     {    return !( p1 == p2 );    }

     private:

     pp_read left_, right_;
     double gap_, dev_;

};

inline void BinaryWrite( int fd, const vec<pp_pair>& h )
{    BinaryWriteComplex( fd, h );    }
inline void BinaryRead( int fd, vec<pp_pair>& h )
{    BinaryReadComplex( fd, h );    }

class pp_mpair : public pp_pair {

     public:

     pp_mpair( const pp_pair& p ) : pp_pair(p), left_mark_(0), 
          right_mark_( p.Right( ).size( ) ) { }

     pp_mpair( const pp_read& left, const pp_read& right, int m1, int m2, 
          double gap, double dev ) : left_mark_(m1), right_mark_(m2)
     {    SetLeft(left);
          SetRight(right);
          SetGap(gap);
          SetDev(dev);    }

     int LeftMark( ) const { return left_mark_; }
     int RightMark( ) const { return right_mark_; }
     void SetLeftMark( int m ) { left_mark_ = m; }
     void SetRightMark( int m ) { right_mark_ = m; }

     friend Bool operator==( const pp_mpair& p1, const pp_mpair& p2 );
     friend Bool operator<( const pp_mpair& p1, const pp_mpair& p2 );

     friend ostream& operator<<( ostream& o, const pp_mpair& p );

     private:

     int left_mark_, right_mark_;

};

typedef istring pp_closure;

class HyperKmerPath; // forward declaration

void ClosePairedPairs( const vec<pp_pair>& pairs, const vec<Bool>& pairs_to_close,
     const vec<int>& L, const double dmult, vec< vec<pp_closure> >& closures, 
     vec< vec<double> >& devs, int max_ext, vec<Bool>& fail, int max_processed = 0, 
     int max_unprocessed = 0, int verbosity = 0, int simple_walk_verbosity = 0,
     int max_opens = 0, int max_nodes = 0, Bool create_closures_if_fail = False,
     const int max_closures = 0 );

void ClosePairedPairsEasy( int min_side, const HyperKmerPath& h,
     const vec<pp_pair>& pairs, const vec<Bool>& pairs_to_close, const vec<int>& L, 
     const double dmult, vec< vec<pp_closure> >& closures, vec< vec<double> >& devs, 
     int max_ext, vec<Bool>& fail, int max_processed, int max_unprocessed, 
     int verbosity, int simple_walk_verbosity, int max_opens = 0, 
     int max_nodes = 0, Bool create_closures_if_fail = False,
     const int max_closures = 0 );

void ExtendPairedPairs( const vec<pp_pair>& pairs, const vec<int>& L, 
     const double dmult, vec<pp_pair>& extensions );

void GetOverlaps( const pp_read& r1, const pp_read& r2, vec<int>& offsets,
     const Bool require_first_subset_second = False );
void GetOverlaps( const pp_read& r1, const pp_read& r2, vec<int>& offsets, 
     const vec<int>& L, int min_overlap );

Bool ReadyToClose( const pp_pair& p, const double dmult );

void GetClosures( const pp_mpair& p, const vec<int>& L, const double dmult,
     vec<pp_closure>& closures, const Bool trim = True, const Bool verbose = False );

void FindJoins( const pp_mpair& p, const pp_pair& x, const vec<int>& L,
     const double dmult, vec<pp_mpair>& pnew, int wt1 = 0, int wt2 = 0,
     Bool must_extend = True );

void JoinReads( const pp_read& r1, const pp_read& r2, const int offset, pp_read& r );

void FindJoins(
     /* inputs: */  const pp_read& r1, const double p1, const double d1,
                    const pp_read& r2, const double p2, const double d2,
                    const double dmult, const vec<int>& L,
     /* outputs: */ vec<pp_read>& r, vec<double>& p, vec<double>& d,
                    vec<int>& offsets );

int OverlapLength( const pp_read& r1, const pp_read& r2, const int offset,     
     const vec<int>& L );

int ClosePairedPairsDirectUnique( vec<pp_pair>& ppp, const HyperKmerPath& h,
     int max_loops, const vec<Bool>& remove );

Bool IsClosed( const pp_pair& p, const vec<int>& L );

Bool ValidOffset( const pp_read& x, const pp_read& y, int o );

#endif
