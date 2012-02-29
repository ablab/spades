///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef BLOCK_H
#define BLOCK_H

#include "Alignment.h"
#include "math/Functions.h"
#include "system/System.h"



/*
 * Class block
 *
 * It captures the information of an aligning block from one align: this
 * is a gapless aligning span, identified uniquely by start on read1,
 * start on read2, and length.
 */
class block {
  
public:
  
  block( ) :
    len_ ( 0 ), mut_ ( 0 ), begin1_ ( 0 ), begin2_ ( 0 ) { }
  
  block( int len, int mut, int begin1, int begin2 ) :
    len_ ( len ), mut_ ( mut ), begin1_ ( begin1 ), begin2_ ( begin2 ) { }
  
  void Set( int len, int mut, int begin1, int begin2 )
    { len_ = len; mut_ = mut; begin1_ = begin1; begin2_ = begin2; }
  
  void SetLen( int len ) { len_ = len; }
  
  void SetMut( int mut ) { mut_ = mut; }
  
  void SetBegin1( int begin1 ) { begin1_ = begin1; }
  
  void SetBegin2( int begin2 ) { begin2_ = begin2; }
  
  int Len( ) const { return len_; }

  int Mut( ) const { return mut_; }
  
  int Begin1( ) const { return begin1_; }
  
  int Begin2( ) const { return begin2_; }
  
  int End1( ) const { return begin1_ + len_; }
  
  int End2( ) const { return begin2_ + len_; }
  
  friend ostream &operator<< ( ostream &out, const block &the_block ) {
    out << the_block.len_ << "\t"
	<< the_block.mut_ << "\t"
	<< the_block.begin1_ << "\t"
	<< the_block.begin2_ << "\n";
    return out; }
  
  friend istream &operator>> ( istream &in, block &the_block ) {
    int len;
    int mut;
    int begin1;
    int begin2;
    in >> len >> mut >> begin1 >> begin2;
    the_block.Set( len, mut, begin1, begin2 );
    return in; }


private:
  
  int len_;            // length of block
  int mut_;            // number of mutations in block
  int begin1_;         // begin of block on read1
  int begin2_;         // begin of block on read2

};



/*
 * BlocksOverlap
 * 
 * Returns the amount of overlap between blocks. An example: say that
 * read1 and read2 are chimp reads, and that they align the same human
 * chromosme as per blocks left and right. This means that the two
 * intervals on human are I1= [left.begin2_, left.begin2_ + len_], and
 * I2 = [right.begin2_, right.begin2_ + len_]. BlocksOverlap returns
 * the overlap between I1 and I2 with a twist: if either block's
 * mutation rate is too high, then pretend the two do not overlap. 
 */
inline int BlocksOverlap( const block &left, const block &right )
{
  const float max_mut_rate = 0.4;
  float left_mut_rate = SafeQuotient( left.Mut( ), left.Len( ) );
  float right_mut_rate = SafeQuotient( right.Mut( ), right.Len( ) );

  if (  left_mut_rate > max_mut_rate || right_mut_rate > max_mut_rate )
    return 0;
  
  int begin = Max( left.Begin2( ), right.Begin2( ) );
  int end = Min( left.End2( ), right.End2( ) );
  
  return Max( 0, end - begin );
}



#endif
