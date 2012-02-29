///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file CharString.cc
 * \author ghall
 */
#include "feudal/CharString.h"
#include "feudal/FeudalStringDefs.h"
#include "feudal/OuterVec.h"
#include "system/ErrNo.h"
#include <cmath>
#include <cstring>

using std::cout;
using std::endl;

/*
 * External Broad Requirements
 */

template class SmallVec< char, MempoolAllocator< char > >;
template class OuterVec< String >;
template class FeudalString< char >;

// ToStringAbbrev
String ToStringAbbrev(const long x)
{ if (abs(x) < 1000) return ToString(x);
  double xDouble = x;
  xDouble /= 1000.;
  String suffix = "K";
  unsigned int precision = 0;
  if (fabs(xDouble) > 1000.)
  { xDouble /= 1000.;
    suffix = "M";
    ++precision; }
  if (fabs(xDouble) > 1000.)
  { xDouble /= 1000.;
    suffix = "G";
    ++precision; }
  return ToString(xDouble, precision) + suffix; }

// ToStringAddCommas
String ToStringAddCommas(const long x)
{ String tmps = x < 0 ? ToString(-x) : ToString(x);
  String result;
  result.reserve(tmps.size()+tmps.size()/3);
  String::const_reverse_iterator end(tmps.crend());
  unsigned int phase = 0;
  for ( String::const_reverse_iterator itr(tmps.crbegin()); itr != end; ++itr )
  {
      if ( phase == 3 ) { result.push_back(','); phase = 0; }
      phase += 1;
      result.push_back(*itr);
  }
  if ( x < 0 ) result.push_back( '-' );
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

// ToLower
String ToLower(const String &s)
{ String slower(s); slower.ToLower(); return slower; }

// ToUpper
String ToUpper(const String & s)
{ String supper(s); supper.ToUpper(); return supper; }

// DeleteLeadingWhiteSpace
void DeleteLeadingWhiteSpace(String& s)
{ size_type lastpos = 0;
  while((lastpos < s.size()) && isspace(s[lastpos]))
    ++lastpos;
  s.erase(0, lastpos); }

// DeleteTrailingWhiteSpace
void DeleteTrailingWhiteSpace(String& s)
{ size_type firstpos = s.size();
  while((firstpos > 0) && isspace(s[firstpos - 1]))
    --firstpos;
  s.erase(firstpos); }

// WhiteSpaceFree
String WhiteSpaceFree(const String& s)
{ String answer;
  answer.reserve(s.size());
  String::const_iterator end(s.end());
  for ( String::const_iterator itr(s.begin()); itr != end; ++itr )
      if ( !isspace(*itr) ) answer.push_back(*itr);
  return answer; }

// cmp_numeric
bool cmp_numeric(const String& s1, const String& s2)
{ for (String::size_type i = 0; i < s1.size( ); i++)
  { if (i == s2.size()) return false;
    if (!isdigit(s1[i]) || !isdigit(s2[i]))
    { if (s1[i] < s2[i]) return true;
      if (s1[i] > s2[i]) return false;
      continue; }
    String::size_type j1, j2;
    for (j1 = i + 1; j1 < s1.size(); j1++)
      if (!isdigit(s1[j1])) break;
    for (j2 = i + 1; j2 < s2.size(); j2++)
      if (!isdigit(s2[j2])) break;
    long n1 = 0, n2 = 0;
    for (String::size_type k = i; k < j1; k++)
      n1 = (10 * n1) + (s1[k] - '0');
    for (String::size_type k = i; k < j2; k++)
      n2 = (10 * n2) + (s2[k] - '0');
    if (n1 < n2) return true;
    if (n1 > n2) return false;
    if (j1 < j2) return false;
    if (j1 > j2) return true;
    i = j1 - 1; }
    return s1.size() < s2.size(); }

// BaseAlpha
String BaseAlpha(unsigned int n)
{ String result;
  result.reserve(7);
  do
  {
      result.push_back('A' + n%26);
  }
  while ( n /= 26 );
  using std::reverse;
  reverse(result.begin(),result.end());
  return result; }

unsigned int UnBaseAlpha(const String& s)
{ unsigned int answer = 0;
  String::const_iterator end(s.end());
  for ( String::const_iterator itr(s.begin()); itr != end; ++itr )
      answer = 26*answer + (*itr - 'A');
  return answer; }

void BinaryWrite(int fd, const String& s)
{ String::size_type len = s.size();
  if (write(fd, &len, sizeof(len)) != sizeof(len))
  { ErrNo err;
    cout << "BinaryWrite of String length failed" << err << endl;
    CRD::exit(1); }
  BinaryWriteContent(fd, s); }

void BinaryWriteContent(int fd, const String& s)
{ String::size_type len = s.size();
  if (write(fd, s.data(), len) != len)
  { ErrNo err;
    cout << "BinaryWrite of String data failed" << err << endl;
    CRD::exit(1); } }

void BinaryRead(int fd, String& s)
{ String::size_type len;
  if (read(fd, &len, sizeof(len)) != sizeof(len))
  { ErrNo err;
    cout << "BinaryRead of String length failed" << err << endl;
    CRD::exit(1); }
  s.resize(len);
  BinaryReadContent(fd, s); }

void BinaryReadContent(int fd, String& s)
{ String::size_type len = s.size();
  if (read(fd, &s[0], len) != len)
  { ErrNo err;
    cout << "BinaryRead of String data failed" << err << endl;
    CRD::exit(1); } }

std::istream& operator>> (std::istream& in, String& s)
{ std::string buf; in >> buf; s = buf.c_str(); return in; }
