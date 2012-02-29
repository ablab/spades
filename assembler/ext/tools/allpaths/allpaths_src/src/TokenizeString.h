///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "String.h"
#include "Vec.h"

#ifndef TOKENIZE_STRING_H
#define TOKENIZE_STRING_H


/// Tokenize a string (separators: " ", and "\t").
int Tokenize( const String &a_string,
	      vec<String> &tokens );

/// Tokenize a string (use one separator).
int Tokenize( const String &a_string,
              const char sep,
	          vec<String> &tokens );

/// Tokenize a string (use given separators).
int Tokenize( const String &a_string,
	      const vec<char> &separators,
	      vec<String> &tokens);

///Split a string into a vector of strings at a delimiter.
///Any container supporting clear() and push_back(String) will work.
template<class VEC>
void Tokenize(String in, VEC & result, const String & delim=" ") {
  result.clear();
  for (int pos = in.Position(delim); pos!=-1; pos = in.Position(delim)) {
    result.push_back(in.substr(0,pos));
    in.erase(0,pos+delim.size());
  }
  if (!in.empty()) {
    result.push_back(in);
  }
}

/// Tokenize a string (use given separators and don't merge them).
/// E.g. "this,is,comma,,separated" -> "this" "is" "comma" "" "separated"
int TokenizeStrictly( const String &a_string,
		      const vec<char> &separators,
		      vec<String> &tokens);

/// Implement some regex matching by forking off to perl.
void RegexMatch(String regex, String &str, vec<String> &matches);

#endif //TOKENIZE_STRING_H
