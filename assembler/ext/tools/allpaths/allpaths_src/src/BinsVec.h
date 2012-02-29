/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_BinsVec_h
#define __INCLUDE_BinsVec_h

#include <algorithm>
#include "Vec.h"
#include "math/Functions.h"

/**
   File: BinsVec.h

   Defines classes that divide a numeric line into bins and associate a value with each bin.

   @file
*/

// Type: bin_id_t
// Identifies a data bin.  For N bin boundaries there are N+1 bins.
typedef size_t bin_id_t;

/**
   Class: BinsVec

   An array that divides the number line of a numeric type into bins according to given
   bin boundaries, and associates a value to each bin.  The bins can be accessed by
   specifying an index value on the number line; the index value is mapped into the appropriate bin,
   and a reference to the data value stored in the bin is returned.  This is done by
   using <function call operator>.
   The data value in the bins can also be accessed by <bin id> using <operator[]>.

   This class inherits from <vec>, but please *do not change the size* of the vector!!
   That is, do not call any methods of vec that might change the size of the vector;
   the size is fixed and is determined by the bin boundaries passed at construction time.

   Template arguments:

      IDX - a numeric type for the values used to index the array.  Must be a model of the
         STL concept LessThanComparable.  The number line of this numeric type is divided
	 into N+1 bins by N bin boundaries specified at construction time; in each bin, _one_
	 value of the VAL type is stored.  The value can be accessed by <bin id> as well as
	 by an index value that is mapped onto a bin.

      VAL - the type of value stored in each bin.  Must be a type that can be stored in an STL
         vector, that is, must be a model of the Assignable STL concept.

   See also: <BinsVec2>, <Histogram>, <PredictionStats>.  See <UnipathEval> for usage example.
*/
template <class IDX, class VAL>
class BinsVec: public vec<VAL> {
 public:
  typedef vec<VAL> PARENT;
  typedef typename PARENT::size_type bin_id_t;

  BinsVec(): firstBinLo_( IDX() ), lastBinHi_( IDX() ), binBoundariesSet_( False ) { }

  /**
     Constructor: BinsVec constructor
    
     Create a BinsVec with the given bin boundaries.
    
     Parameters:
    
        binBoundaries - the bin boundaries.  *Must be sorted in ascending order!*
          See <ParseIntSet()> and <ParseDoubleSet()>.
  */
  BinsVec(const vec<IDX>& binBoundaries):
    PARENT(binBoundaries.size()+1), binBoundaries_(binBoundaries),
    firstBinLo_( IDX() ), lastBinHi_( IDX() ), binBoundariesSet_( False ) {
    Assert( isValid() );
  }


  /**
     Constructor: BinsVec constructor with default value
    
     Create a BinsVec with the given bin boundaries, and fill the data in each bin
     with the specified value.  
    
     Parameters:
    
        binBoundaries - the bin boundaries.  *Must be sorted in ascending order!*
          See <ParseIntSet()> and <ParseDoubleSet()>.
        defaultVal - the default value to put into each bin.
  */
  BinsVec(const vec<IDX>& binBoundaries, const VAL& defaultVal):
    PARENT(binBoundaries.size()+1, defaultVal), binBoundaries_(binBoundaries),
    firstBinLo_( IDX() ), lastBinHi_( IDX() ), binBoundariesSet_( False ) {
    Assert( isValid() );
  }

  /// Operator: function call operator
  /// Given an index value, find the appropriate bin for that value and return a reference
  /// to the data value stored in that bin.
  typename PARENT::reference operator()( const IDX& i ) { return (*this)[FindBin(i)]; }
  typename PARENT::const_reference operator()( const IDX& i ) const { return (*this)[FindBin(i)]; }

  /// Method: GetNumBins
  /// Returns the number of bins.
  bin_id_t GetNumBins() { return PARENT::size(); }

  /// Method: GetBinLo
  /// Returns the lower bound of a bin, given the <bin id>.  The lower bound of the lowest bin is
  /// *undefined* until it is set by <SetOuterBounds>.  (The bin receives all values below the
  /// lowest <bin boundary>).
  const IDX& GetBinLo(bin_id_t binId) const {
    Assert( binBoundariesSet_ ) ;
    return binId == 0 ? firstBinLo_ : binBoundaries_[binId-1];
  }

  /// Method: GetBinHi
  /// Returns the uppder bound of a bin, given the <bin id>.  The upper bound of the highest bin is
  /// *undefined* until it is set by <SetOuterBounds>.  (The bin receives all values above and
  /// including the highest <bin boundary>).
  const IDX& GetBinHi(bin_id_t binId) const {
    Assert( binBoundariesSet_ ) ;
    return binId == binBoundaries_.size()  ? lastBinHi_ : binBoundaries_[binId];
  }

  /// Method: FindBin
  /// Given an index value, find the appropriate bin for it and return the <bin id>.
  /// Note that you can directly get a reference to the contents of the appropriate bin
  /// by using the <function call operator>.  Note also that the <bin boundaries>
  /// passed at construction time *must* be sorted in ascending order and duplicate-free
  /// for this to work.
  bin_id_t FindBin( const IDX& idx ) const {
    return Min( (bin_id_t)
		distance( binBoundaries_.begin(),
			  lower_bound( binBoundaries_.begin(), binBoundaries_.end(), idx ) ),
		(bin_id_t) binBoundaries_.size() );
  }

  /// Method: SetOuterBounds
  /// Record the lower bound of the lowest bin, and the upper bound of the
  /// highest bin.  The lowest bin receives all values below the lowest <bin boundary>,
  /// and the highest bin receives all values above the highest bin boundary.
  /// This method must be called before calling <GetBinLo()> or <GetBinHi()>.
  void SetOuterBounds(const IDX& firstBinLo, const IDX& lastBinHi) {
    Assert( firstBinLo < binBoundaries_[0]  &&  lastBinHi > binBoundaries_[binBoundaries_.size()-1] );
    firstBinLo_ = firstBinLo; lastBinHi_ = lastBinHi;  binBoundariesSet_ = True;
  }
  
 private:
  /// Private field: binBoundaries
  /// The boundaries of the bins.  For n boundaries there are n+1 bins.
  /// The boundaries denote half-open intervals: a bin includes the value
  /// at its lower boundary but not at its upper boundary.  The lowest bin
  /// has no lower boundary and just receives all values below the lowest bin
  /// boundary; for display purposes, call <SetOuterBounds()> to specify
  /// the lower boundary of the lowest bin.  Analogously for the highest bin.
  vec<IDX> binBoundaries_;

  /// Private field: firstBinLo_
  /// The lower boundary of the lowest bin.  Used for display purposes only;
  /// the bin simply receives all values below the lowest <bin boundary>.
  /// Set by <SetOuterBounds()>.
  IDX firstBinLo_;

  /// Private field: lastBinHi_
  /// The upper bin boundary.  Used for display purposes only;
  /// the bin simply receives all values above and including the highest <bin boundary>.
  /// Set by <SetOuterBounds()>.
  IDX lastBinHi_;

  /// Private field: binBoundariesSet_
  /// Whether <SetOuterBounds()> has been called to set the outer boundaries of the bins.
  Bool binBoundariesSet_;

 protected:

  /// Protected method: isValid
  /// Tests whether the representation of this object is valid / internally consistent.
  /// Used to test the validity of constructor arguments.
  Bool isValid() const {
    ForceAssert( binBoundaries_.nonempty() );
    ForceAssertEq( binBoundaries_.size()+1, PARENT::size() );
    ForceAssert( is_sorted( binBoundaries_.begin(), binBoundaries_.end() ) );
      // there are no duplicate bin boundaries
    ForceAssert( adjacent_find( binBoundaries_.begin(), binBoundaries_.end() ) == binBoundaries_.end() );
    return True;
  }
  
};  // class BinsVec

/**
   Class: BinsVec2

   An array that divides a number plane (two orthogonal number lines of numeric types)
   into square bins according to given
   bin boundaries, and associates a value to each bin.  The bins can be accessed by
   specifying an index value on each number line; the pair of index values is mapped into the appropriate bin,
   and a reference to the data value stored in the bin is returned.  This is done by
   using <function call operator>.
   The data value in the bins can also be accessed by supplying a pair of <bin ids>,
   giving the coordinate of the bin on each of the two axes, using <operator[]>.

   This class inherits from <vec>, but please *do not change the size* of the vector!!
   That is, do not call any methods of vec that might change the size of the vector;
   the size is fixed and is determined by the bin boundaries passed at construction time.

   Template arguments:

      IDX1, IDX2 - a pair of numeric types for the values used to index the array.  Each must be a model of the
         STL concept LessThanComparable.  The number line of each numeric type is divided
	 into bins by bin boundaries specified at construction time; in each square bin (corresponding to a pair
	 of linear bins, one on each number line), _one_
	 value of the VAL type is stored.  The value can be accessed by a pair of <bin ids> as well as
	 by a pair of index values that is mapped onto a bin.

      VAL - the type of value stored in each bin.  Must be a type that can be stored in an STL
         vector, that is, must be a model of the Assignable STL concept.

   See also: <BinsVec>, <Histogram>, <PredictionStats>.  See <UnipathEval> for usage example.
*/
template <class IDX1, class IDX2, class VAL>
class BinsVec2: public BinsVec< IDX1, BinsVec<IDX2, VAL> > {
 public:
  typedef BinsVec< IDX1, BinsVec<IDX2, VAL> > PARENT;
  typedef typename PARENT::bin_id_t bin_id_t;  

  /**
     Constructor: BinsVec2 constructor
    
     Create a BinsVec with the given bin boundaries.
    
     Parameters:
    
        binBoundaries1, binBoundaries2 - the bin boundaries on the two
	   number lines.  Each *must be sorted in ascending order, with no duplicates!*
	   See <ParseIntSet()> and <ParseDoubleSet()>.
  */
  BinsVec2(const vec<IDX1>& binBoundaries1, const vec<IDX2>& binBoundaries2):
    PARENT(binBoundaries1, BinsVec<IDX2, VAL>( binBoundaries2 ) ) {  }

  /**
     Constructor: BinsVec2 constructor with default value
    
     Create a BinsVec2 with the given bin boundaries, and fill the data in each bin
     with the specified value.  
    
     Parameters:
    
        binBoundaries1, binBoundaries2 - the bin boundaries on the two number lines.
	*Must be sorted in ascending order!*
	See <ParseIntSet()> and <ParseDoubleSet()>.
	
        defaultVal - the default value to put into each bin.
  */
  BinsVec2(const vec<IDX1>& binBoundaries1, const vec<IDX2>& binBoundaries2, const VAL& defaultVal):
    PARENT(binBoundaries1, BinsVec<IDX2, VAL>( binBoundaries2, defaultVal ) ) { }

  typename BinsVec<IDX2, VAL>::reference operator()( const IDX1& i1, const IDX2& i2 ) { return (*(PARENT *)this)(i1)(i2); }
  typename BinsVec<IDX2, VAL>::const_reference operator()( const IDX1& i1, const IDX2& i2 ) const { return ((*(PARENT *)this)(i1))(i2); }

  /// Method: SetOuterBounds
  /// Record the lower bound of the lowest bin, and the upper bound of the
  /// highest bin, on each of the two number lines.
  /// The lowest bin receives all values below the lowest <bin boundary>,
  /// and the highest bin receives all values above the highest bin boundary.
  /// This method must be called before calling <GetBinBounds()>.
  void SetOuterBounds(const IDX1& firstBinLo1, const IDX1& lastBinHi1, const IDX2& firstBinLo2, const IDX2& lastBinHi2) {
    PARENT::SetOuterBounds( firstBinLo1, lastBinHi1 );
    for (typename PARENT::size_type i=0; i<PARENT::size(); i++)
      (*this)[i].SetOuterBounds( firstBinLo2, lastBinHi2 );
  }
  

  /// Method: GetNumBins2
  /// Return the number of bins on the second number axis.
  bin_id_t GetNumBins2() { retun (*this)(0).GetNumBins(); }

  /// Method: GetBinBounds
  /// Return the bounds of the bin identified by a pair of <bin ids>.
  void GetBinBounds(bin_id_t binId1, bin_id_t binId2,
			  IDX1& binLo1, IDX1& binHi1, IDX2& binLo2, IDX2& binHi2) const {
    const BinsVec<IDX2, VAL>& binsVec2 = (*(PARENT *)this)(0);
    binLo1 = PARENT::GetBinLo(binId1);
    binHi1 = PARENT::GetBinHi(binId1);
    binLo2 = binsVec2.GetBinLo(binId2);
    binHi2 = binsVec2.GetBinHi(binId2);
  }
  
};  // class BinsVec2



#endif
// #ifndef __INCLUDE_BinsVec_h

// Synonyms: Various synonyms
//
//    bin id - See <bin_id_t>
//    bin boundary - See <binBoundaries_>

