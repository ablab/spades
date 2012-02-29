///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////


//  Id: Mutmer.h,v 1.4 2000/07/12 16:30:30 jaffe Exp $             

#ifndef MUTMER
#define MUTMER

// ================================================================================
//
// The mutmer class stores four integers: pos1, pos2, len, and errors.  
//
// Everything that follows relates to the INTERPRETATION of mutmers, and has no
// reflection in the class definition....
//
// A mutmer is supposed to represent a string of length len from the alphabet 
// {A,C,G,T,N}, which matches two reads at positions pos1 and pos2.  The
// integer errors counts the number of N's, which should correspond to the positions
// in the match where the two reads differ.
//
// We also assume that a mutmer satisfies the following two conditions:
// (a) N does not appear on either end;
// (b) There are at least two bases between any two N's.
//
// ================================================================================

#include "Basevector.h"

class mutmer {

     public:

     mutmer( ) { }

     mutmer(int pos1, int pos2, int len, int e) :
          pos1_(pos1), pos2_(pos2), len_(len), e_(e) { }

     void SetFrom(int pos1, int pos2, int len, int e)
     {    pos1_ = pos1;
          pos2_ = pos2;
          len_ = len;
          e_ = e;    }
               
     void SetFrom( const mutmer &orig)
     {    pos1_ = orig.pos1_;
          pos2_ = orig.pos2_;
          len_ = orig.len_;
          e_ = orig.e_;    }
               
     void Unpack(int& pos1, int& pos2, int& len, int& e) const
     {    pos1 = pos1_;
          pos2 = pos2_;
          len = len_;
          e = e_;    }

     int Pos1() const { return pos1_; }
     int Pos2() const { return pos2_; }
     int End1() const { return pos1_ + len_; }
     int End2() const { return pos2_ + len_; }
     int Offset() const { return pos1_ - pos2_; }
     int Length() const { return len_; }
     int Errors() const { return e_; }

     void SetPos1( const int pos1 ) { pos1_ = pos1; }
     void SetPos2( const int pos2 ) { pos2_ = pos2; }
     void SetLength( const int len ) { len_ = len; }
     void SetErrors( const int errs ) { e_ = errs; }
     
     void CountErrors( const basevector &rd1, const basevector &rd2 )
       {
	 e_ = 0;
	 for ( int i = 0; i < len_; ++i )
	   if ( rd1[ pos1_ + i ] != rd2[ pos2_ + i ] )
	     ++e_;
       }

     void TrimLeft( const int bases, const basevector &rd1, const basevector &rd2 )
       {
	 // cout << "Pre left trim by " << bases << ": " << endl;
	 // this->Print( cout, rd1, rd2 );
	 if ( bases >= len_ )
	 {
	   len_ = 0;
	   e_ = 0;
	 }
	 else
	 {
	   for ( int i = 0; i < bases; ++i )
	     if ( rd1[ pos1_ + i ] != rd2[ pos2_ + i ] )
	       --e_;
	   AssertGe( e_, 0 );
	   pos1_ += bases;
	   pos2_ += bases;
	   len_ -= bases;
	 }
	 // cout << "Post left trim by " << bases << ": " << endl;
	 // this->Print( cout, rd1, rd2 );
       }	   

     void TrimRight( const int bases, const basevector &rd1, const basevector &rd2 )
       {
	 // cout << "Pre right trim by " << bases << ": " << endl;
	 // this->Print( cout, rd1, rd2 );
	 if ( bases >= len_ )
	 {
	   len_ = 0;
	   e_ = 0;
	 }
	 else
	 {
	   len_ -= bases;
	   for ( int i = 0; i < bases; ++i )
	     if ( rd1[ pos1_ + len_ + i ] != rd2[ pos2_ + len_ + i ] )
	       --e_;
	   AssertGe( e_, 0 );
	 }
	 // cout << "Post right trim by " << bases << ": " << endl;
	 // this->Print( cout, rd1, rd2 );
       }	   

     friend bool operator< ( const mutmer &lhs, const mutmer &rhs )
       {
	 if ( lhs.pos1_ < rhs.pos1_ ) return true;
	 if ( lhs.pos1_ > rhs.pos1_ ) return false;
	 if ( lhs.pos2_ < rhs.pos2_ ) return true;
	 if ( lhs.pos2_ > rhs.pos2_ ) return false;
	 if ( lhs.len_ < rhs.len_ ) return true;
	 if ( lhs.len_ > rhs.len_ ) return false;
	 if ( lhs.e_ < rhs.e_ ) return true;
	 return false;
       }

     // Print a summary of the mutmer.
     friend ostream & operator<< ( ostream &out, const mutmer &mm )
       {
	 return out << mm.Pos1() << "-" << mm.End1() << "\t"
		    << mm.Pos2() << "-" << mm.End2() << "\t"
		    << "(" << mm.Offset() << ")\t"
		    << mm.Length() << "\t"
		    << mm.Errors();
       }

     // Print a human-readable version of the mutmer with the aligning
     // bases from each read.
     void Print( ostream &out, const basevector &rd1, const basevector &rd2 );

     private:

     int pos1_, pos2_, len_, e_;

};

#include <functional>

struct mutmer_is_short : public unary_function<mutmer,bool>
{
  int min_length;

  mutmer_is_short( int min ) { min_length = min; }

  bool operator() ( const mutmer &m ) const
    { return m.Length() < min_length; }
};

#endif
