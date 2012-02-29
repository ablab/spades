///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PARSE_RANGE_H
#define PARSE_RANGE_H

#include "String.h"

template<typename T>
T ConvertFromString(const String& value);


template<typename T>
bool ParseRange(const String range, T& lower, T& upper,
		bool& lower_inc, bool& upper_inc,
		bool& lower_open, bool& upper_open) {

  // Need a range string at least 3 charactes long with a comma.
  if (range.empty())
    return false;
  if (range.Freq(",") != 1)
    return false;
  if (range.size() < 3)
    return false;

  // Extract range start and stop modifiers
  String lowerChar = String(range[0]);
  String upperChar = String(range[range.size() - 1]);

  // Check if lower value in range is inclusive or exclusive
  if (lowerChar == "[")
    lower_inc = true;
  else if (lowerChar == "(" || lowerChar == "]")
    lower_inc = false;
  else
    return false;

  // Check if upper value in range is inclusive or exclusive
  if (upperChar == "]")
    upper_inc = true;
  else if (upperChar == ")" || upperChar == "[")
    upper_inc = false;
  else
    return false;

  // Extract upper and lower range values as strings
  String lowerStr = range.Between(lowerChar,",");
  String upperStr = range.Between(",", upperChar);

  // Check if range is open, otherwise extract lower value
  if (lowerStr.empty())
    lower_open = true;
  else {
    lower_open = false;
    lower = ConvertFromString<T>(lowerStr);
  }

  // Check if range is open, otherwise extract upper value
  if (upperStr.empty())
    upper_open = true;
  else {
    upper_open = false;
    upper = ConvertFromString<T>(upperStr);
  }

  // We have a valid range!
  return true;
}

#endif
