///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SEQ_INTERVAL_H
#define SEQ_INTERVAL_H

#include <functional>
#include <algorithm>
#include <iostream>

using namespace std;

/*
 * class seq_interval
 *
 * A sequence interval class. It is a (simply connected) interval in
 * a given but generic genomic sequence. E.g.: a read in a contig, a
 * contig in a supercontig, an insert in a supercontig.
  // TODO: Potentially dangerous truncation of IDs
 */
class seq_interval {

public:

  seq_interval( );

  seq_interval( int interval_id, int seq_id, int begin, int end );

  void Set( int interval_id, int seq_id, int begin, int end );

  void SetIntervalId( int interval_id ) { interval_id_ = interval_id; }

  void SetSeqId( int seq_id ) { seq_id_ = seq_id; }

  void SetBegin( int begin ) { begin_ = begin; }

  void SetEnd( int end ) { end_ = end; }

  void IncrementEnd() { ++end_; }

  int IntervalId( ) const { return interval_id_; }

  int SeqId( ) const { return seq_id_; }

  int Begin( ) const { return begin_; }

  int End( ) const { return end_; }

  int Length( ) const { return end_ - begin_; }

  // Set the coordinates of this to be the intersection of a and b.
  // If a and b are on different sequences, this is set to be the
  // interval of zero length beginning and ending at the beginning of
  // a.  In all cases, the interval id is untouched.
  void SetToIntersectionOf( const seq_interval &a, const seq_interval &b );

  // If this and other have different seq_id_ return false, otherwise
  // return true and fill gap with the gap size.
  bool GapSizeWith( const seq_interval &other, int &gap ) const;

  // Does this interval overlap the other interval?
  bool HasOverlapWith( const seq_interval &other ) const;

  // What is the amount of overlap with the other interval?
  int HasAmountOfOverlapWith( const seq_interval &other ) const;

  // operator==
  friend bool operator== ( const seq_interval &ll, const seq_interval &rr ) {
    if ( ll.interval_id_ != rr.interval_id_ ) return false;
    if ( ll.seq_id_ != rr.seq_id_ ) return false;
    if ( ll.begin_ != rr.begin_ ) return false;
    if ( ll.end_ != rr.end_ ) return false;
    return true;
  }

  // operator!=
  friend bool operator!= (const seq_interval &lint, const seq_interval &rint) {
    return ( ! ( lint == rint ) );
  }

  // Sort first by sequence id, then start on sequence.
  friend bool operator< (const seq_interval &lint, const seq_interval &rint) {
    return ( lint.seq_id_ < rint.seq_id_ ||
	     ( lint.seq_id_ == rint.seq_id_ &&
	       ( lint.begin_ < rint.begin_ ||
		 ( lint.begin_ == rint.begin_ &&
		   ( lint.end_ < rint.end_ ||
		     ( lint.end_ == rint.end_ &&
		       ( lint.interval_id_ < rint.interval_id_ ) ) ) ) ) ) );
  }

  // operator>
  friend bool operator> (const seq_interval &lint, const seq_interval &rint) {
    if ( lint < rint ) return false;
    if ( lint == rint ) return false;
    return true;
  }

  // Sort by interval_id first.
  struct OrderByIntervalId :
    public binary_function<seq_interval,seq_interval,bool> {
    bool operator()(const seq_interval& lint, const seq_interval& rint) {
      return ( lint.interval_id_ < rint.interval_id_ ||
	       ( lint.interval_id_ == rint.interval_id_ &&
		 ( lint.seq_id_ < rint.seq_id_ ||
		   ( lint.seq_id_ == rint.seq_id_ &&
		     ( lint.begin_ < rint.begin_ ||
		       ( lint.begin_ == rint.begin_ &&
			 ( lint.end_ < rint.end_ ) ) ) ) ) ) );
    }
  };

  // Sort by length.  In the event of a tie, sort by interval_id.
  struct OrderByLength :
    public binary_function<seq_interval,seq_interval,bool> {
    bool operator()(const seq_interval& lint, const seq_interval& rint) {
      return ( lint.Length() < rint.Length() ||
	       ( lint.Length() == rint.Length() &&
		 ( lint.interval_id_ < rint.interval_id_ ||
		   ( lint.interval_id_ == rint.interval_id_ &&
		     ( lint.seq_id_ < rint.seq_id_ ||
		       ( lint.seq_id_ == rint.seq_id_ &&
			 ( lint.begin_ < rint.begin_ ) ) ) ) ) ) );
    }
  };

  friend ostream& operator<< ( ostream &out, const seq_interval &seq_int );

  friend istream& operator>> ( istream &in, seq_interval &seq_int );


private:

  // TODO: potentially dangerous truncation of indexes
  int interval_id_; // id of interval
  int seq_id_;      // id of the sequence containing the interval
  int begin_;       // begin base on sequence
  int end_;         // end base on sequence

};

#endif
