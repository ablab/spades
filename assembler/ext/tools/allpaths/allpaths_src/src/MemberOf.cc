///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Vec.h"
#include "ParseSet.h"
#include "MemberOf.h"


bool IsMemberOfSet(const String& set, const String& value) {
  vec<String> validValues;
  ParseStringSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const int value) {
  vec<int> validValues;
  ParseIntSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const longlong value) {
  vec<longlong> validValues;
  ParseLongLongSet(set, validValues, false);
  return Member(validValues, value);
}

bool IsMemberOfSet(const String& set, const unsigned int value) {
  int valueLongLong = static_cast<longlong> (value);
  return IsMemberOfSet(set, valueLongLong);
}

bool IsMemberOfSet(const String& set, const double value) {
  vec<double> validValues;
  ParseDoubleSet(set, validValues, false);
  return Member(validValues, value);
}
