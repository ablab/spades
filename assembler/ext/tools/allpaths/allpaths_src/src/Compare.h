///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file Compare.h
 * \author tsharpe
 * \date Mar 4, 2010
 *
 * \brief
 */
#ifndef COMPARE_H_
#define COMPARE_H_

#include <algorithm>
#include <string>
#include <vector>

template <class T>
int compare( T const& t1, T const& t2 )
{
    return (t2 < t1) - (t1 < t2);
}

template <class T, class A>
int compare( std::vector<T,A> const& v1, std::vector<T,A> const& v2 )
{
    typedef typename std::vector<T,A>::const_iterator Itr;
    Itr itr1 = v1.begin();
    Itr itr2 = v2.begin();
    using std::min; Itr end = itr1 + min(v1.size(),v2.size());
    int result = 0;
    while ( !result && itr1 != end )
    { result = ::compare(*itr1,*itr2); ++itr1; ++itr2; }
    if ( !result ) result = ::compare(v1.size(),v2.size());
    return result;
}

inline int compare( std::string const& s1, std::string const& s2 )
{ return s1.compare(s2); }

template <class T>
struct Comparator
{
    int operator()(T const& t1, T const& t2)
    { return compare(t1,t2); }
};

#endif /* COMPARE_H_ */
