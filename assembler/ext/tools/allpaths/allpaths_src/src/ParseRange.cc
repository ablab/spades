///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "ParseRange.h"

template<> int ConvertFromString<int>(const String& value) {
  return value.Int();
}

template<> unsigned int ConvertFromString<unsigned int>(const String& value) {
  return value.Int();
}

template<> long ConvertFromString<long>(const String& value) {
  return value.Int(); // Int() returns a long
}

template<> double ConvertFromString<double>(const String& value) {
  return value.Double();
}

template<> String ConvertFromString<String>(const String& value) {
  return value;
}
