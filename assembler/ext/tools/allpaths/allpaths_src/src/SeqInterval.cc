///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "SeqInterval.h"
#include "VecTemplate.h"

BINARY2_DEF(seq_interval);



/*
 * seq_interval
 * Constructor
 */
seq_interval::seq_interval( ) :
  interval_id_ ( -1 ),
  seq_id_ ( -1 ),
  begin_ ( -1 ),
  end_ ( -1 )
{ }



/*
 * seq_interval
 * Constructor
 */
seq_interval::seq_interval( int interval_id, int seq_id, int begin, int end ) :
  interval_id_ ( interval_id ),
  seq_id_ ( seq_id ),
  begin_ ( begin ),
  end_ ( end )
{ }

/*
 * seq_interval
 * Set
 */
void seq_interval::Set( int interval_id, int seq_id, int begin, int end )
{
  interval_id_ = interval_id;
  seq_id_ = seq_id;
  begin_ = begin;
  end_ = end;
}

/*
 * seq_interval
 * SetToIntersectionOf
 */
void seq_interval::SetToIntersectionOf( const seq_interval &a, const seq_interval &b )
{
  if ( a.seq_id_ != b.seq_id_ )
  {
    seq_id_ = a.seq_id_;
    begin_ = a.begin_;
    end_ = a.begin_;
  }

  seq_id_ = a.seq_id_;
  begin_ = max( a.begin_, b.begin_ );
  end_ = min( a.end_, b.end_ );

  if ( end_ < begin_ )
    end_ = begin_;
}

/*
 * seq_interval
 * GapSizeWith
 */
bool seq_interval::GapSizeWith( const seq_interval &other, int &gap ) const
{
  gap = other.Begin( ) - this->End( );
  return ( this->SeqId( ) == other.SeqId( ) );
}

/*
 * seq_interval
 * HasOverlapWith
 */
bool seq_interval::HasOverlapWith( const seq_interval &other ) const
{
  if ( HasAmountOfOverlapWith( other ) )
    return true;
  else
    return false;
}

/*
 * seq_interval
 * HasAmountOfOverlapWith
 */
int seq_interval::HasAmountOfOverlapWith( const seq_interval &other ) const
{
   if ( seq_id_ != other.SeqId())
     return 0;

   int begin = max(begin_, other.Begin());
   int end = min(end_, other.End());
   return max(0, end - begin);
}


/*
 * seq_interval
 * operator<<
 */
ostream& operator<< ( ostream &out, const seq_interval &seq_int )
{
  out << seq_int.interval_id_ << "\t"
      << seq_int.seq_id_ << "\t"
      << seq_int.begin_ << "\t"
      << seq_int.end_;
  
  return out;
}



/*
 * seq_interval
 * operator>>
 */
istream& operator>> ( istream &in, seq_interval &seq_int )
{
  in >> seq_int.interval_id_
     >> seq_int.seq_id_
     >> seq_int.begin_
     >> seq_int.end_;
  
  return in;
}



