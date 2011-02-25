#include "toyexamples.h"

#include <string>
#include <set>
#include <map>
#include <iostream>
#include <cassert>
#include <fstream>

//#include "../graphVisualizer.hpp"

using namespace std;
//using namespace gvis;

typedef map < string, unsigned > CKMerSet;


class CSubstringIterator
{
  string const s;
  unsigned const substr_size;
  unsigned cur_pos;
  bool Last;

public:
  CSubstringIterator ( string const & s, unsigned substr_size ) :
    s ( s ),
    substr_size ( substr_size ),
    cur_pos ( 0 ),
    Last ( false )
  {
    assert ( s.size () >= substr_size );
  }

  void Reset   () { cur_pos = 0; Last = false; }
  void Advance () { ++cur_pos; if ( s.size () == cur_pos ) Last = true; }
  bool IsLast  () const { return Last; }

  string operator * () const
  {
    if ( cur_pos + substr_size <= s.size () )
      return s.substr ( cur_pos, substr_size );

    string prefix ( s.begin () + cur_pos, s.end () );
    return prefix += s.substr ( 0, substr_size - prefix.size() );
  }
};



extern void ConstructDeBruijnGraph ( string genome, unsigned k )//, string filename )
{
  CKMerSet S;

  // preparing a file
  ofstream f;
  f.open ( "graph.dot" );
  f << "digraph " << genome << " {" << endl;

  // adding vertices
  CSubstringIterator it ( genome, k );
  for ( ; ! it.IsLast (); it.Advance () )
  {
    if ( ( S.insert ( make_pair ( *it, S.size () ) ) ).second )
      f << ( S.size () - 1 ) << " [label" << "=" << *it << "]" << endl;
  }

  // adding edges
  CSubstringIterator it1 ( genome, k + 1 );
  for ( ; ! it1.IsLast (); it1.Advance () )
  {
    string const m = *it1;
    CKMerSet::const_iterator from = S.find ( m.substr ( 0, k ) );
    CKMerSet::const_iterator to   = S.find ( m.substr ( 1, k ) );
    assert ( from != S.end () && to != S.end () );
    f << from -> second << "->" << to -> second << endl;

    //[penwidth=5,color=red]
  }


  f << "}" << endl;
  f.close ();

  return;
}
