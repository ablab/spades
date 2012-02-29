///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Class alignlet used in scaffolding.

#ifndef ALIGNLET_H
#define ALIGNLET_H

#include "lookup/LookAlign.h"

// An alignlet is a compact version of a look_align.

class alignlet {

public:

  alignlet( ) { }

  alignlet( const look_align& la ) {
    pos2_ = la.pos2( ), Pos2_ = la.Pos2( );
    target_id = la.target_id;
    target_length_ = la.Fw1( ) ? la.target_length : - la.target_length - 1;
  }

  alignlet( const int pos2, const int Pos2, const int target_id,
       const int target_length, const Bool fw1 ) 
       : pos2_(pos2), Pos2_(Pos2), target_id(target_id)
  {    target_length_ = fw1 ? target_length : -target_length - 1;    }
  
  int pos2( ) const { return pos2_; }
  int Pos2( ) const { return Pos2_; }
  int Fw1( ) const { return target_length_ >= 0; }
  void Shift(int offset) { pos2_ += offset; Pos2_ += offset; }
  int TargetId( ) const { return target_id; }
  int TargetLength( ) const {
    return target_length_ < 0 ? - target_length_ - 1 : target_length_;
  }

  void SetTargetId( const int id ) { target_id = id; }
  void SetTargetLength( const int length, const bool fw ) {
    target_length_ = ( fw ? length : - length - 1 );
  }

  void Reverse( ) {
    target_length_ = - target_length_ - 1;
    int xpos2 = this->TargetLength( ) - Pos2( );
    int xPos2 = this->TargetLength( ) - pos2( );
    pos2_ = xpos2;
    Pos2_ = xPos2;
  }

  void PrintReadableBrief( ostream& out, const String& query_name ) const {
    out << query_name << ( target_length_ < 0 ? "rc" : "fw" ) << " vs " << target_id << " from " << pos2_ << " to " << Pos2_ 
	<< " (of " << ( target_length_ < 0 ? - target_length_ -1 : target_length_ ) << ")\n";

  }
  
private:

  int pos2_, Pos2_;
  int target_id;
  int target_length_;    // signed: encapsulates orientation (fw1)

};

#endif
