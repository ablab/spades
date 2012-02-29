/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// VecOverlap provides a solution to the following problem: given vec<T> objects
// x1,...,xn, build a database so that if given vec<T> objects y1,...,ym,
// one can in turn look each yi up in S, returning all proper overlaps between
// yi and an xj.

#ifndef VEC_OVERLAP_H
#define VEC_OVERLAP_H

#include "CoreTools.h"

template<class T> class vec_overlap {

 public:
  
  // The constructor builds the data structures needed to compute overlaps.
  
  vec_overlap( const vec< vec<T> >& x );

  ~vec_overlap();
  
  // GetOverlaps returns a sorted list of pairs (j,o) where j is an index to one 
  // of the x's and o is the offset in the alignment, e.g.
  // ---------y---------
  //  ----------xi----------
  // has an offset of 1.
  //
  // If allow_subset_x = False, then xi has to extend past y on at
  // least one end.
  
  void GetOverlaps( const vec<T>& y, vec< pair<int,int> >& overlaps,
                    Bool allow_subset_x = True ) const;
  
 private:
  class imp;

  imp* m_pImp;
};

#endif
