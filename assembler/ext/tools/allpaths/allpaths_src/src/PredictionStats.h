/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef __INCLUDE_PredictionStats_h
#define __INCLUDE_PredictionStats_h

#include <iostream>
#include "String.h"

const bool PREDICTED_NO = false;
const bool PREDICTED_YES = true;
const bool ACTUAL_NO = false;
const bool ACTUAL_YES = true;

/**
   Class: PredictionStats

   A 4x4 table showing statistics about our success in predicting
   some binary property of objects in a set of objects (for example,
   which kmers are <genomic kmers> or which unipaths are <unique unipaths>).
*/
class PredictionStats {
 public:
  PredictionStats(const String& propertyYes = "Yes", const String& propertyNo = "No"):
    propertyYes_(propertyYes), propertyNo_(propertyNo) {
    (*this)( PREDICTED_NO, ACTUAL_NO ) =
      (*this)( PREDICTED_NO, ACTUAL_YES ) =
      (*this)( PREDICTED_YES, ACTUAL_NO ) =
      (*this)( PREDICTED_YES, ACTUAL_YES ) = 0;
  }

  const size_t& operator()(bool predicted, bool actual) const
    { return counts_[bool2int(predicted)][bool2int(actual)]; }
  size_t& operator() (bool predicted, bool actual)
    { return counts_[bool2int(predicted)][bool2int(actual)]; }

  const String& PropertyYesName() const { return propertyYes_; }
  const String& PropertyNoName() const { return propertyNo_; }

  void Reset( const String& propertyYes, const String& propertyNo ) {
    propertyYes_ = propertyYes;
    propertyNo_ = propertyNo;
    (*this)( PREDICTED_NO, ACTUAL_NO ) =
      (*this)( PREDICTED_NO, ACTUAL_YES ) =
      (*this)( PREDICTED_YES, ACTUAL_NO ) =
      (*this)( PREDICTED_YES, ACTUAL_YES ) = 0;
  }

  void Print( ostream& s = cout ) const {
    s << "   predicted " << propertyYes_ << " actual " << propertyYes_ << ": " << (*this)( PREDICTED_YES, ACTUAL_YES )  << endl;
    s << "   predicted " << propertyNo_ << " actual " << propertyYes_ << ": " << (*this)( PREDICTED_NO, ACTUAL_YES )  << endl;
    s << "   predicted " << propertyYes_ << " actual " << propertyNo_ << ": " << (*this)( PREDICTED_YES, ACTUAL_NO )  << endl;
    s << "   predicted " << propertyNo_ << " actual " << propertyNo_ << ": " << (*this)( PREDICTED_NO, ACTUAL_NO )  << endl;
  }

 private:
  size_t counts_[2][2];
  String propertyYes_, propertyNo_;

  static inline int bool2int(bool b) { return b ? 1 : 0; }
};

inline ostream& operator<< ( ostream& s, const PredictionStats& ps ) {
  ps.Print( s );
  return s;
}

#endif
// #ifndef __INCLUDE_PredictionStats_h

