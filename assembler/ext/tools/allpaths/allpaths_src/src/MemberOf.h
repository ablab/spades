///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MEMBER_OF_H
#define MEMBER_OF_H

#include "String.h"
#include "ParseSet.h"
#include "ParseRange.h"

template<typename T>
bool IsMemberOfRange(const T lower, const T upper, const T& value,
		     const bool lower_inc=true, const bool upper_inc=true,
		     const bool lower_open=false, const bool upper_open= false) {
  bool valid = true;
  valid &= lower_open
    || (lower_inc && value >= lower)
    || (!lower_inc && value > lower);
  valid &= upper_open
    || (upper_inc && value <= upper)
    || (!upper_inc && value < upper);
  return valid;
}

template<typename T>
bool IsMemberOfRange(const String& range, const T value) {
  T lower, upper;
  bool lower_inc, upper_inc;
  bool lower_open, upper_open;
  if (ParseRange(range, lower, upper, lower_inc, upper_inc,
		 lower_open, upper_open))
    return IsMemberOfRange(lower, upper, value, lower_inc, upper_inc,
			   lower_open, upper_open);
  else
    FatalErr("Invalid Range");
}

template<typename T>
bool IsMemberOf(const String& valid, const T value) {
  if (!valid.empty()) {
    if (valid[0] == '(' || valid[0] == '[' || valid[0] == ']')
      return IsMemberOfRange(valid, value);
    else
      return IsMemberOfSet(valid, value);
  }
  return false;
}


bool IsMemberOfSet(const String& valid, const String& value);

bool IsMemberOfSet(const String& valid, const int value);

bool IsMemberOfSet(const String& valid, const unsigned int value);

bool IsMemberOfSet(const String& valid, const longlong value);

bool IsMemberOfSet(const String& valid, const double value);


#endif

