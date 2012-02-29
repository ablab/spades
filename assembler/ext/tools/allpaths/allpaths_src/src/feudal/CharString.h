///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file CharString.h
 * \author ghall
 * \date Oct 1, 2009
 *
 * \brief Instantiation of FeudalString as char
 */
#ifndef __CHARSTRING_H__
#define __CHARSTRING_H__

#include "feudal/FeudalString.h"
#include "math/Calculator.h"
#include <iomanip> /* for setprecision */
#include <sstream>

typedef FeudalString<char, std::char_traits<char> > String;

// ToString type conversion suite
inline String ToString(int x)
{ char cstr[12]; sprintf(cstr, "%d", x); return String(cstr); }
inline String ToString(unsigned int x)
{ char cstr[12]; sprintf(cstr, "%u", x); return String(cstr); }
inline String ToString(long x)
{ char cstr[21]; sprintf(cstr, "%ld", x); return String(cstr); }
inline String ToString(unsigned long x)
{ char cstr[21]; sprintf(cstr, "%lu", x); return String(cstr); }
inline String ToString(float x, int precision)
{ std::ostringstream oss; oss.setf(std::ios::fixed);
  oss << std::setprecision(precision) << x; return oss.str(); }
inline String ToString(double x, int precision)
{ std::ostringstream oss; oss.setf(std::ios::fixed);
  oss << std::setprecision(precision) << x; return oss.str(); }
inline String ToString(char const* cstr) { return cstr; }
inline String ToString(const String& s) { return s; }
inline String ToString(const std::string& s) { return s.c_str(); }
template <class T>
inline String ToString(const T& sumthin)
{ std::ostringstream oss; oss << sumthin; return oss.str(); }

// ToStringBool
inline String ToStringBool(bool value) { return (value ? "True" : "False"); }
inline String ToString(bool value) { return ToStringBool(value); }

// Format a possibly large integer for human consumption by putting it in terms
// of thousands, millions, or billions, with appropriate suffix.
String ToStringAbbrev(const long x);

// Format a possibly large integer for human consumption by putting commas
// between each 3-digit group.
String ToStringAddCommas(const long x);

// Transform a string to lowercase using isupper and tolower.  Not safe for
// international strings, which is probably irrelevant.
String ToLower(const String & s);

// Transform a string to lowercase using islower and toupper.  Not safe for
// international strings, which is probably irrelevant.
String ToUpper(const String & s);

// Remove leading white space from a string
void DeleteLeadingWhiteSpace(String& s);

// Remove trailing white space from a string
void DeleteTrailingWhiteSpace(String& s);

// Remove leading and trailing white space from a string
String WhiteSpaceFree(const String& s);

// cmp_numeric: compare two strings, in such a way that they will be sorted
// numerically if they contain integers.  For example, chr1 < chr2 < chr10.
bool cmp_numeric(const String& s1, const String& s2);

// BaseAlpha: convert a nonnegative integer into an alphabetic string,
// according to the order A, B, ..., Z, AA, AB, ..., etc.
String BaseAlpha(unsigned int n);
unsigned int UnBaseAlpha(const String& s);

// Feudal File Access
void BinaryRead(int fd, String& s);
void BinaryWrite(int fd, const String& s);
void BinaryReadContent(int fd, String& s);
void BinaryWriteContent(int fd, const String& s);

// Evaluate a mathematical expression
inline String Evaluate(const String& expr, const bool /* verbose */ = 0)
{ return ToString(evaluate(expr.c_str())); }

inline String operator+ (const std::string& lhs, const String& rhs)
{ return lhs.c_str() + rhs; }
inline String operator+ (const String& lhs, const std::string& rhs)
{ return lhs + rhs.c_str(); }

inline bool operator== (const std::string& lhs, const String& rhs)
{ return lhs.c_str() == rhs; }
inline bool operator== (const String& lhs, const std::string& rhs)
{ return lhs == rhs.c_str(); }

inline bool operator!= (const std::string& lhs, const String& rhs)
{ return lhs.c_str() != rhs; }
inline bool operator!= (const String& lhs, const std::string& rhs)
{ return lhs != rhs.c_str(); }

inline bool operator> (const std::string& lhs, const String& rhs)
{ return lhs.c_str() > rhs; }
inline bool operator> (const String& lhs, const std::string& rhs)
{ return lhs > rhs.c_str(); }

inline bool operator>= (const std::string& lhs, const String& rhs)
{ return lhs.c_str() >= rhs; }
inline bool operator>= (const String& lhs, const std::string& rhs)
{ return lhs >= rhs.c_str(); }

inline bool operator< (const std::string& lhs, const String& rhs)
{ return lhs.c_str() < rhs; }
inline bool operator< (const String& lhs, const std::string& rhs)
{ return lhs < rhs.c_str(); }

inline bool operator<= (const std::string& lhs, const String& rhs)
{ return lhs.c_str() <= rhs; }
inline bool operator<= (const String& lhs, const std::string& rhs)
{ return lhs <= rhs.c_str(); }

std::istream& operator>> (std::istream& in, String& s);
inline std::ostream& operator<< (std::ostream& out, const String& s)
{ return out << s.c_str(); }

/// stream operators

/// similar to getline for std::string, but treats '\r' '\n' and "\r\n" as
/// delimiters
inline std::istream& getline(std::istream& in, String& s)
{ String tmp;
  char c;
  // TODO: like the previous implementation, this is flawed in that it leaves
  // the istream in a failed state when returning a final line that lacks a
  // delimiter.
  while ( in.get(c) )
  { if ( c == '\r' && in.peek() == '\n' ) in.get(c);
    if ( c == '\n' || c == '\r' ) break;
    tmp.push_back(c); }
  s = tmp;
  return in; }

/// return a string filled by reading to the end of a stream
inline String Slurp(std::istream& is)
{ String tmp; char c; while (is.get(c)) tmp.push_back(c); return tmp; }

#endif // __CHARSTRING_H__
