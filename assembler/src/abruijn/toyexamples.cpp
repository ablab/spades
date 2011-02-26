#include "toyexamples.hpp"
#include "hash.hpp"

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



extern void ConstructDeBruijnGraph ( string genome, unsigned read_size, unsigned k, )//, string filename )
{
  CKMerSet S;

  // preparing a file
  ofstream f;
  f.open ( "graph.dot" );
  f << "## genome=" << genome << " readsize=" << read_size << " k=" << k << endl;
  f << "digraph " << genome << " {" << endl;

  //////////////////////////// DE BRUIJN GRAPH
  // adding vertices
  for ( CSubstringIterator kmer_it ( genome, k ); ! kmer_it.IsLast (); kmer_it.Advance () )
  {
    if ( ( S.insert ( make_pair ( *kmer_it, S.size () ) ) ).second )
      f << ( S.size () - 1 ) << " [label" << "=" << *kmer_it << "]" << endl;
  }

  // adding edges
  for ( CSubstringIterator edge_it ( genome, k + 1 ); ! edge_it.IsLast (); edge_it.Advance () )
  {
    string const edge = *edge_it;
    CKMerSet::const_iterator from = S.find ( edge.substr ( 0, k ) );
    CKMerSet::const_iterator to   = S.find ( edge.substr ( 1, k ) );
    assert ( from != S.end () && to != S.end () );
    f << from -> second << "->" << to -> second << endl;
  }



  ///////////////////////////// A BRUIJN GRAPH
  for ( CSubstringIterator read_it ( genome, read_size ); ! read_it.IsLast (); read_it.Advance () )
  {
    string const read = *read_it;
    HashSym < string > h;
    multimap < unsigned, string > kmers_with_hash;
    for ( CSubstringIterator kmer_it ( read, k ); ! kmer_it.IsLast (); kmer_it.Advance () )
      kmers_with_hash.insert ( make_pair ( h ( *kmer_it ), *kmer_it ) );
    multimap < unsigned, string >::const_iterator tmp_it = kmers_with_hash.begin ();
    CKMerSet::const_iterator from = S.find ( tmp_it -> second );
    assert ( from != S.end () )
    while ( tmp_it -> second == kmers_with_hash.begin () -> second )
      ++tmp_it;
    if ( tmp_it != kmers_with_hash.end () ) 
    // this may happen if all the k-mers of the current read are equal
    {
      CKMerSet::const_iterator to = S.find ( tmp_it -> second );
      f << from -> second << "->" << to -> second << " [penwidth=5,color=red]" << endl;
    }
  }


  f << "}" << endl;
  f.close ();

  return;
}
