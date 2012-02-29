///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

/*
 * \file FeudalString.h
 * \author ghall
 * \date Oct 28, 2009
 *
 * \brief Definitions for FeudalString
 */

#ifndef FEUDAL_STRING_DEF_H_
#define FEUDAL_STRING_DEF_H_

#include "feudal/FeudalString.h"
#include "feudal/SmallVecDefs.h"
#include "feudal/OuterVecDefs.h"
#include "system/Assert.h"
#include "system/System.h"
#include <cctype>
#include <cstddef>
#include <cstring>
#include <sstream>

template <class charT, class Traits>
const charT FeudalString<charT, Traits>::EmptyString[] = "";


// c_str()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::const_pointer
FeudalString<charT, Traits>::c_str() const
{ if (empty()) return EmptyString;
  container &vec = const_cast<container&>(mContainer);
  vec.reserve(size() + 1);
  *vec.dataEnd() = charT();
  return data(); }

// find()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ AssertLe(pos,size());
  if ( n <= size() )
  { const_pointer end = mContainer.dataEnd() - n;
    for ( const_pointer itr = data()+pos; itr <= end; ++itr )
      if ( !Traits::compare(itr,cstr,n) )
        return itr - data(); }
  return npos; }


// find()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find(FeudalString<charT, Traits>::value_type c,
                              FeudalString<charT, Traits>::size_type pos) const
{ AssertLe(pos,size());
  const_pointer ppp = Traits::find(data()+pos,size()-pos,c);
  return ppp ? ppp-data() : npos; }


// rfind()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::rfind(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ if ( pos == npos ) pos = size();
  AssertLe(pos,size());
  if ( n <= size() )
  { using std::min;
    const_pointer end = data();
    const_pointer itr = end + min(pos,size()-n) + 1;
    while ( itr > end )
      if ( !Traits::compare(--itr,cstr,n) )
        return itr - end; }
  return npos; }

// rfind()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::rfind(
        FeudalString<charT, Traits>::value_type c,
        FeudalString<charT, Traits>::size_type pos) const
{ if ( pos == npos ) pos = size();
  AssertLe(pos,size());
  if ( size() )
  { using std::min; pos = min(pos,size()-1);
    const_pointer end = data();
    const_pointer itr = end + min(pos,size()-1) + 1;
    while ( itr > end )
      if ( Traits::eq(*--itr,c) )
        return itr - end; }
  return npos; }

// find_first_of()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_first_of(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ AssertLe(pos,size());
  const_pointer end = mContainer.dataEnd();
  for ( const_pointer itr = data()+pos; itr < end; ++itr )
      if ( Traits::find(cstr,n,*itr) )
          return itr - data();
  return npos; }

// find_last_of
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_last_of(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ if ( pos == npos ) pos = size();
  AssertLe(pos,size());
  if ( size() )
  { const_pointer end = data();
    const_pointer itr = end + min(pos,size()-1) + 1;
    while ( itr > end )
      if ( Traits::find(cstr,n,*--itr) )
        return itr - end; }
  return npos; }

// find_first_not_of
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_first_not_of(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ AssertLe(pos,size());
  const_pointer end = mContainer.dataEnd();
  for ( const_pointer itr = data()+pos; itr < end; ++itr )
      if ( !Traits::find(cstr,n,*itr) )
          return itr - data();
  return npos; }

// find_first_not_of
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_first_not_of(
        typename FeudalString<charT, Traits>::value_type c,
        typename FeudalString<charT, Traits>::size_type pos) const
{ AssertLe(pos,size());
  const_pointer end = mContainer.dataEnd();
  for ( const_pointer itr = data()+pos; itr < end; ++itr )
      if ( !Traits::eq(c,*itr) )
          return itr - data();
  return npos; }


// find_last_not_of
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_last_not_of(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos,
        typename FeudalString<charT, Traits>::size_type n) const
{ if ( pos == npos ) pos = size();
  AssertLe(pos,size());
  if ( size() )
  { const_pointer end = data();
    const_pointer itr = end + min(pos,size()-1) + 1;
    while ( itr > end )
      if ( !Traits::find(cstr,n,*--itr) )
        return itr - end; }
  return npos; }

// find_last_not_of
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::find_last_not_of(
        typename FeudalString<charT, Traits>::value_type c,
        typename FeudalString<charT, Traits>::size_type pos) const
{ if ( pos == npos ) pos = size();
  AssertLe(pos,size());
  if ( size() )
  { const_pointer end = data();
    const_pointer itr = end + min(pos,size()-1) + 1;
    while ( itr > end )
      if ( !Traits::eq(*--itr,c) )
        return itr - end; }
  return npos; }

// compare()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::compare(
        typename FeudalString<charT, Traits>::const_pointer cstr) const
{ size_t mySize = size();
  size_t cstrSize = Traits::length(cstr);
  using std::min;
  int result = Traits::compare(data(), cstr, min(mySize,cstrSize));
  if ( !result ) result = ::compare(mySize,cstrSize);
  return result; }

// compare()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::compare(
        typename FeudalString<charT, Traits>::size_type pos1,
        typename FeudalString<charT, Traits>::size_type n1,
        typename FeudalString<charT, Traits>::const_pointer cstr) const
{ size_t mySize = lenLimit(pos1,n1);
  size_t cstrSize = Traits::length(cstr);
  using std::min;
  int result = Traits::compare(data()+pos1, cstr, min(mySize,cstrSize));
  if ( !result ) result = ::compare(mySize,cstrSize);
  return result; }

// compare()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::compare(
        typename FeudalString<charT, Traits>::size_type pos1,
        typename FeudalString<charT, Traits>::size_type n1,
        const FeudalString<charT, Traits>& str) const
{ n1 = lenLimit(pos1,n1);
  using std::min;
  int result = Traits::compare(data()+pos1, str.data(), min(n1, str.size()));
  if ( !result ) result = ::compare(n1,str.size());
  return result; }

// compare()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::compare(
        typename FeudalString<charT, Traits>::size_type pos1,
        typename FeudalString<charT, Traits>::size_type n1,
        const FeudalString<charT, Traits>& str,
        typename FeudalString<charT, Traits>::size_type pos2,
        typename FeudalString<charT, Traits>::size_type n2) const
{ n1 = lenLimit(pos1,n1);
  n2 = str.lenLimit(pos2,n2);
  using std::min;
  int result = Traits::compare(data()+pos1, str.data()+pos2, min(n1, n2));
  if ( !result ) result = ::compare(n1,n2);
  return result; }

// compare()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::compare(
        typename FeudalString<charT, Traits>::size_type pos1,
        typename FeudalString<charT, Traits>::size_type n1,
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type n2) const
{ n1 = lenLimit(pos1,n1);
  using std::min;
  int result = Traits::compare(data()+pos1, cstr, min(n1, n2));
  if ( !result ) result = ::compare(n1,n2);
  return result; }

// IsInt()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::IsInt( long* pVal ) const
{ std::istringstream iss(String(c_str(),size()));
  long val;
  if ( !(iss >> val) ) return false; // can't extract long
  char c1, c2;
  if ( (iss >> c1) )
  { if ( (iss >> c2) ) return false; // two trailing non-numeric characters
    if ( c1 == 'K' ) val *= 1000L;
    else if ( c1 == 'M' ) val *= 1000000L;
    else if ( c1 == 'G' ) val *= 1000000000L;
    else if ( c1 == '%' ) val /= 100L;
    else return false; } // can't interpret trailing char as multiplier
  if ( pVal ) *pVal = val;
  return true; }

// Int()
template <class charT, class Traits>
long
FeudalString<charT, Traits>::Int() const
{ long val;
  if (!IsInt(&val))
    FatalErr("Attempt to convert string \"" + *this + "\" to an int failed.");
  return val; }

// IsDouble()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::IsDouble( double* pVal ) const
{ std::istringstream iss(String(c_str(),size()));
  double val;
  if ( !(iss >> val) ) return false; // can't extract double
  char c1, c2;
  if ( (iss >> c1) )
  { if ( (iss >> c2) ) return false; // two trailing non-double characters
    if ( c1 == 'K' ) val *= 1000.;
    else if ( c1 == 'M' ) val *= 1000000.;
    else if ( c1 == 'G' ) val *= 1000000000.;
    else if ( c1 == '%' ) val /= 100.;
    else return false; } // can't interpret trailing char as multiplier
  if ( pVal ) *pVal = val;
  return true; }

// Double()
template <class charT, class Traits>
double
FeudalString<charT, Traits>::Double() const
{ double val;
  if (!IsDouble(&val))
    FatalErr("Attempt to convert string \""<<*this<<"\" to a double failed.");
  return val; }

// IsBool()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::IsBool( bool* pVal ) const
{ char const* cstr = c_str();
  if ( !cstr[0] || !memcmp(cstr,"False",6) )
  { if ( pVal ) *pVal = false; return true; }
  if ( !memcmp(cstr,"True",5) )
  { if ( pVal ) *pVal = true; return true; }
  return false; }

// ToBool()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::ToBool() const
{ bool val;
  if (!IsBool(&val) )
    FatalErr("Attempt to convert string \"" << *this << "\" to a bool failed.");
  return val; }

// Toupper()
template <class charT, class Traits>
FeudalString<charT, Traits>&
FeudalString<charT, Traits>::ToUpper()
{ iterator stop(end());
  for ( iterator itr(begin()); itr != stop; ++itr )
      *itr = toupper(*itr);
  return *this; }

// ToLower()
template <class charT, class Traits>
FeudalString<charT, Traits>&
FeudalString<charT, Traits>::ToLower()
{ iterator stop(end());
  for ( iterator itr(begin()); itr != stop; ++itr )
      *itr = tolower(*itr);
  return *this; }

// Position()
template <class charT, class Traits>
int
FeudalString<charT, Traits>::Position(const FeudalString<charT, Traits>& str,
        typename FeudalString<charT,Traits>::size_type endSearchAt) const
{ int_size_assert();
  AssertLe(endSearchAt,size());
  using std::min; endSearchAt = min(endSearchAt,size());
  size_type srchLen = str.size();
  if ( endSearchAt >= srchLen )
  { const_pointer cstr = str.data();
    const_pointer end = data() + endSearchAt - srchLen;
    for ( const_pointer itr = data(); itr <= end; ++itr )
      if ( !Traits::compare(itr,cstr,srchLen) )
        return itr - data(); }
  return npos; }

// Contains()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::Contains(
        typename FeudalString<charT, Traits>::const_pointer cstr,
        typename FeudalString<charT, Traits>::size_type pos) const
{ size_type len = Traits::length(cstr);
  using std::min; if ( pos == npos ) pos = size() - min(size(),len);
  return lenLimit(pos,len) == len && !Traits::compare(data()+pos,cstr,len); }

// Contains()
template <class charT, class Traits>
bool
FeudalString<charT, Traits>::Contains(const FeudalString<charT, Traits>& str,
        typename FeudalString<charT, Traits>::size_type pos) const
{ size_type len = str.size();
  using std::min; if ( pos == npos ) pos = size() - min(size(),len);
  return lenLimit(pos,len) == len &&
          !Traits::compare(data()+pos,str.data(),len); }

// Before()
template <class charT, class Traits>
FeudalString<charT, Traits>
FeudalString<charT, Traits>::Before(
        typename FeudalString<charT, Traits>::const_pointer cstr) const
{ size_type idx = find(cstr);
  if (idx == npos)
    FatalErr("Before(): Couldn't find \""<<cstr<<"\" in \""<< *this << "\".\n");
  return substr(0, idx); }

// Before()
template <class charT, class Traits>
FeudalString<charT, Traits>
FeudalString<charT, Traits>::Before(
        const FeudalString<charT, Traits>& str) const
{ size_type idx = find(str);
  if (idx == npos)
    FatalErr("Before(): Couldn't find \"" + str + "\" in \"" + *this + "\".\n");
  return substr(0, idx); }

// After()
template <class charT, class Traits>
FeudalString<charT, Traits>
FeudalString<charT, Traits>::After(
        typename FeudalString<charT, Traits>::const_pointer cstr) const
{ size_type idx = find(cstr);
  if (idx == npos)
    FatalErr("After() couldn't find \"" <<cstr<< "\" in \"" << *this << "\".");
  return substr(idx + Traits::length(cstr)); }

// After()
template <class charT, class Traits>
FeudalString<charT, Traits>
FeudalString<charT, Traits>::After(const FeudalString<charT, Traits>& str) const
{ size_type idx = find(str);
  if (idx == npos)
    FatalErr("After() couldn't find \"" << str << "\" in \"" << *this << "\".");
  return substr(idx + str.size()); }

// ReplaceBy()
template <class charT, class Traits>
FeudalString<charT, Traits>&
FeudalString<charT, Traits>::ReplaceBy(const FeudalString<charT, Traits>& from,
                                        const FeudalString<charT, Traits>& to)
{ size_type idx = find(from);
  if (idx == npos)
    FatalErr("ReplaceBy: could not find \""<<from<<"\" in \""<< *this << "\".");
  return replace(begin(idx),begin(idx+from.size()),to.begin(),to.end()); }

// GlobalReplaceBy()
template <class charT, class Traits>
FeudalString<charT, Traits>&
FeudalString<charT, Traits>::GlobalReplaceBy(
        const FeudalString<charT, Traits>& from,
        const FeudalString<charT, Traits>& to)
{ size_type fromLen = from.size();
  const_iterator toBegin(to.begin());
  const_iterator toEnd(to.end());
  size_type toLen = to.size();
  size_type idx = 0;
  size_type newIdx;
  while ( (newIdx = find(from,idx)) != npos )
  { replace(begin(newIdx),begin(newIdx+fromLen),toBegin,toEnd);
    idx = newIdx + toLen; }
  return *this; }

// Freq()
template <class charT, class Traits>
typename FeudalString<charT, Traits>::size_type
FeudalString<charT, Traits>::Freq(const FeudalString<charT, Traits>& str) const
{ size_type result = 0;
  size_type idx = ~0u;
  while ( (idx = find(str,++idx)) != npos ) result += 1;
  return result; }

// ReplaceExtension()
template <class charT, class Traits>
FeudalString<charT, Traits>
FeudalString<charT, Traits>::ReplaceExtension(FeudalString<charT, Traits> const& ext, FeudalString<charT, Traits> const& newext) const
{ size_type extLen = ext.size();
  size_type len = size();
  const_pointer me = data();
  if ( len >= extLen && !Traits::compare(me+len-extLen,ext.data(),extLen) )
      len -= extLen;
  return String().reserve(len+newext.size()).assign(me,me+len).append(newext); }

template <class charT, class Traits>
void FeudalString<charT,Traits>::tooBigForInts() const
{ FatalErr("Requesting String's size as an integer, but the size is " << size() << " and it won't fit."); }

#endif
