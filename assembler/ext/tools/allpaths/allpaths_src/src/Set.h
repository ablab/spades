///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// Set.h defines some utilities often used with sets, e.g. Member()
// and i/o operators.

#ifndef SET_H
#define SET_H

#include <iostream>
#include <set>
#include "CoreTools.h"

using namespace std;

template<class T> bool Member( const set<T>& the_set, const T& value )
{    return the_set.find(value) != the_set.end( );    }  // Breaks cxx 

template<class T> ostream& operator<<(ostream& out, const set<T>& the_set)
{    out << the_set.size( ) << "\n";
     typename set<T>::const_iterator the_set_iter = the_set.begin(); 
     for ( ; the_set_iter != the_set.end(); ++the_set_iter )
       out << *the_set_iter << "\n";
     return out;    }

template<class T> istream& operator>>(istream& in, set<T>& s)
{    int n;
     in >> n;
     char c;
     in.get(c);
 
     T temp;
     for ( int idx = 0; idx < n; ++idx )
     {
       in >> temp;
       s.insert( temp );
     }
     return in;    }


/// Get all values in a set
template<class Value, typename Cmp> vec<Value> values(const std::set<Value,Cmp> &  m) {
  vec<Value> ret(m.size());
  int i = 0;
  for (typename std::set<Value,Cmp>::const_iterator iter = m.begin();
       iter != m.end();
       ++i, ++iter) {
    ret[i] = *iter;
  }
  return ret;
}

#endif
