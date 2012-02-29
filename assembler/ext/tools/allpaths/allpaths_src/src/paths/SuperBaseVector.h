// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef SUPERBASEVECTOR
#define SUPERBASEVECTOR

#include "Vec.h"
#include "Basevector.h"
#include "math/Functions.h"

#include <numeric>
#include <functional>

/// A simple class which holds a series of basevectors with
/// (possibly negative) gaps between them.  This is what a 
/// KmerPath logically maps to in sequence space.  The function
/// KmerBaseBroker::ToSequence() returns an object of this type.
///
/// A SuperBaseVector of size() n has n basevectors, accessed as
/// Seq(0)...Seq(n-1), and n-1 gpas, accessed as Gap(0)...Gap(n-2).

class SuperBaseVector {
public:
  // using constructor/destructor defaults

  // Data accessors.  If you want a non-const basevector, copy it yourself.
  const basevector& Seq( int i ) const { return seqs[i]; }
  const basevector& operator[]( int i ) const { return seqs[i]; }
  
  pair<int,int> Gap( int i ) const { return make_pair(mingaps[i],maxgaps[i]); }
  int MinGap( int i ) const { return mingaps[i]; }
  int MaxGap( int i ) const { return maxgaps[i]; }
  int MinGap() const { return accumulate( mingaps.begin(), mingaps.end(), 0 ); }
  int MaxGap() const { return accumulate( maxgaps.begin(), maxgaps.end(), 0 ); }

  int size() const { return seqs.size(); }

  // Replace *this by its reverse complement.
  void ReverseComplement( ) {
    mingaps.ReverseMe();
    maxgaps.ReverseMe();
    seqs.ReverseMe();
    for_each( seqs.begin(), seqs.end(), 
//            This also works -- the ambiguity is in the template resolution
// 	      mem_fun_ref_t<void,basevector>( &basevector::ReverseComplement) );
	      mem_fun_ref( (basevector & (basevector::*)()) 
			   &basevector::ReverseComplement) );
  }

  

  int ReducedLength( ) const
  {    int total = 0;
       for ( int i = 0; i < seqs.isize( ); i++ )
            total += seqs[i].size( );
       return total;    }


  /// You need a pair of ints to specify a base in a SuperBaseVector.
  /// If sbv is a SuperBaseVector and loc is a SuperBaseVector::Loc,
  /// then the base it refers to is sbv[loc.seq][loc.base]
  struct Loc {
    int seq;
    int base;
    Loc() : seq(-1), base(-1) { } // uninitialized
    Loc(int s, int b) : seq(s), base(b) { }
  };


  // Changing the contents:
  // To modify one of these, you must alternately push back seqs and gaps.
  // Feel free to add methods to make this less fascist if you need them.
  void PushSeq( const basevector& seq ) {
    ForceAssertEq( seqs.size(), mingaps.size() );
    seqs.push_back( seq );
  }
  
  void PushGap( pair<int,int> gap ) {
    ForceAssertEq( seqs.size(), mingaps.size()+1 );
    mingaps.push_back( gap.first );
    maxgaps.push_back( gap.second );
  }
  void PushGap( int a, int b ) {
    ForceAssertEq( seqs.size(), mingaps.size()+1 );
    mingaps.push_back( a );
    maxgaps.push_back( b );
  }


  void Print( ostream& out ) {
    if( size()==0 ) {
      out << "(empty SuperBaseVector)";
    }
    Seq(0).Print( out );
    for(int i=1; i<size(); i++) {
      out << "(gap " << MinGap(i-1) << " - " << MaxGap(i-1) << ")\n";
      Seq(i).Print( out );
    }
  }

  // PrintN: print, showing Ns for gaps.

  void PrintN( ostream& out )
  {    int count = 0;
       for ( int i = 0; i < size( ); i++ )
       {    for ( int j = 0; j < (int) Seq(i).size( ); j++ )
            {    if ( count > 0 && count % 80 == 0 ) out << "\n";
                 out << as_base( Seq(i)[j] );
                 ++count;    }
            if ( i < size( ) - 1 )
            {    int g = Max( 0, ( MinGap(i) + MaxGap(i) ) / 2 );
                 for ( int j = 0; j < g; j++ )
                 {    if ( count > 0 && count % 80 == 0 ) out << "\n";
                      out << "N";
                      ++count;    }    }    }
      out << "\n";    }

private:
  vec<basevector> seqs;
  vec<int> mingaps;
  vec<int> maxgaps;
public:
  String name;  // Just for the user's convenience

};





#endif
