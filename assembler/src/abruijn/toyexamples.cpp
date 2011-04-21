#include "toyexamples.hpp"
#include "hash.hpp"

#include <string>
#include <set>
#include <map>
#include <iostream>
#include <cassert>
#include <fstream>

using namespace std;
using namespace hashing;

typedef map < string, unsigned > CKMerSet;
typedef map < pair < string, string >, unsigned > CWeightedEdgeSet;

long special_lexicographic3 ( string const & s )
{
  assert ( s.size () == 3 );
  if ( s == "TAA" ) return 1;
  if ( s == "ACG" ) return 2;
  if ( s == "CGA" ) return 3;
  if ( s == "GAA" ) return 4;
  if ( s == "CTA" ) return 5;
  if ( s == "ACT" ) return 6;
  if ( s == "AAA" ) return 7;
  if ( s == "AAC" ) return 8;

  assert ( false );

  return 10;
}

long standard_lexicographic3 ( string const & s )
{
  //cout << endl << s;
  assert ( s.length () == 3 );
  if ( s == "AAA" ) return 1;
  if ( s == "AAC" ) return 2;
  if ( s == "ACG" ) return 3;
  if ( s == "ACT" ) return 4;
  if ( s == "CGA" ) return 5;
  if ( s == "CTA" ) return 6;
  if ( s == "GAA" ) return 7;
  if ( s == "TAA" ) return 8;

  assert ( false );

  return 10;
}

Hash < string > h;

long my_hash ( string const & s )
{
  //return h ( s );
  return standard_lexicographic3 ( s );
}



class CSubstringIterator
{
  string const s;
  unsigned const substr_size;
  unsigned cur_pos;
  bool const cycle_string;
  bool Last;

public:
  CSubstringIterator ( string const & s, unsigned substr_size, bool cycle_string ) :
    s ( s ),
    substr_size ( substr_size ),
    cur_pos ( 0 ),
    cycle_string ( cycle_string ),
    Last ( false )
  {
    assert ( s.size () >= substr_size );
  }

  void Reset   ()
  {
    cur_pos = 0;
    Last = false;
  }

  void Advance () {
    ++cur_pos;
    if ( cycle_string && s.size () == cur_pos )
      Last = true;
    if ( ! cycle_string && cur_pos + substr_size - 1 == s.size () )
      Last = true;
  }

  bool IsLast  () const { return Last; }

  string operator * () const
  {
    assert ( ! Last );

    if ( cur_pos + substr_size <= s.size () )
      return s.substr ( cur_pos, substr_size );

    assert ( cycle_string );

    string prefix ( s.begin () + cur_pos, s.end () );
    return prefix += s.substr ( 0, substr_size - prefix.size() );
  }
};



extern void ConstructDeBruijnGraph ( string genome, unsigned read_size, unsigned k )//, string filename )
{
  //cout << "\n" << genome;
  //for ( CSubstringIterator kmer_it ( genome, k, false ); ! kmer_it.IsLast (); kmer_it.Advance () )
  //  cout << "\n" << *kmer_it;
  //
  //return;

  CKMerSet S;
  Hash < string > h;

  // preparing a file
  ofstream f;
  //f.open ( "graph.dot" );
  f.open ( ( genome + ".dot" ).c_str () );
  f << "## genome=" << genome << " readsize=" << read_size << " k=" << k << endl;
  f << "##dot -Tjpg " << genome << ".dot -o " << genome << ".jpg" << endl;
  f << "digraph " << genome << " {" << endl;

  //////////////////////////// DE BRUIJN GRAPH
  // adding vertices
  for ( CSubstringIterator kmer_it ( genome, k, true ); ! kmer_it.IsLast (); kmer_it.Advance () )
  {
    if ( ( S.insert ( make_pair ( *kmer_it, S.size () ) ) ).second )
      f << ( S.size () - 1 ) << " [label" << "=" << *kmer_it << "-" << h ( *kmer_it ) << "]" << endl;
  }

  // adding edges
  for ( CSubstringIterator edge_it ( genome, k + 1, true ); ! edge_it.IsLast (); edge_it.Advance () )
  {
    string const edge = *edge_it;
    


    CKMerSet::const_iterator from = S.find ( edge.substr ( 0, k ) );
    CKMerSet::const_iterator to   = S.find ( edge.substr ( 1, k ) );
    assert ( from != S.end () && to != S.end () );
    f << from -> second << "->" << to -> second << endl;
  }



  ///////////////////////////// A BRUIJN GRAPH
  for ( CSubstringIterator read_it ( genome, read_size, true ); ! read_it.IsLast (); read_it.Advance () )
  {
    string const read = *read_it;
    //cout << "\ncur read: " << read << flush;
    multimap < unsigned, string > kmers_with_hash;
    for ( CSubstringIterator kmer_it ( read, k, false ); ! kmer_it.IsLast (); kmer_it.Advance () )
      kmers_with_hash.insert ( make_pair ( h ( *kmer_it ), *kmer_it ) );
    multimap < unsigned, string >::const_iterator tmp_it = kmers_with_hash.begin ();

    //cout << "\n the min k-mer: " << tmp_it -> second << flush;

    CKMerSet::const_iterator from = S.find ( tmp_it -> second );
    assert ( from != S.end () );
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


extern void ConstructDeBruijnGraphSimplified ( string genome, unsigned read_size, unsigned k )
{
  //Hash < string > h;

  ofstream f;
  f.open ( ( genome + ".dot" ).c_str () );
  f << "## genome=" << genome << " readsize=" << read_size << " k=" << k << endl;
  f << "##dot -Tjpg " << genome << ".dot -o " << genome << ".jpg" << endl;
  f << "digraph " << genome << " { rankdir=\"LR\"; node[shape=\"box\"] " << endl;

  CWeightedEdgeSet E;

  ///////////////////////////// DE BRUIJN GRAPH
  for ( CSubstringIterator edge_it ( genome, k + 1, true ); ! edge_it.IsLast (); edge_it.Advance () )
  {
    string const edge = *edge_it;
    string from = edge.substr ( 0, k );
    string to = edge.substr ( 1, k );

    pair < pair < string, string >, unsigned > e = make_pair ( make_pair ( from, to ), 1 );
    pair < CWeightedEdgeSet::iterator, bool > res = E.insert ( e );
    if ( ! res.second )
      ++ (( res.first ) -> second );
  }

  for ( CWeightedEdgeSet::const_iterator eit = E.begin (); eit != E.end (); ++ eit )
    f << eit -> first . first << "->" << eit -> first . second << "[color=grey,label=\"" << eit -> second << "\"]" << endl;    


  ///////////////////////////// A BRUIJN GRAPH
  E.clear ();

  for ( CSubstringIterator read_it ( genome, read_size, true ); ! read_it.IsLast (); read_it.Advance () )
  {
    string const read = *read_it;
    multimap < unsigned, string > kmers_with_hash;
    for ( CSubstringIterator kmer_it ( read, k, false ); ! kmer_it.IsLast (); kmer_it.Advance () )
      kmers_with_hash.insert ( make_pair ( my_hash ( *kmer_it ), *kmer_it ) );
    multimap < unsigned, string >::const_iterator from = kmers_with_hash.begin ();
    multimap < unsigned, string >::const_iterator to = kmers_with_hash.begin ();
    ++to;

    string from_s = from -> second;
    string to_s = to -> second;

    for ( CSubstringIterator kmer_it ( read, k, false ); ! kmer_it.IsLast () && *kmer_it != from_s; kmer_it.Advance () )
      if ( *kmer_it == to_s )
      {
        swap ( from_s, to_s );
         break;
      }


    pair < pair < string, string >, unsigned > e = make_pair ( make_pair ( from_s, to_s ), 1 );
    pair < CWeightedEdgeSet::iterator, bool > res = E.insert ( e );
    if ( ! res.second )
      ++ (( res.first ) -> second );
  }

  for ( CWeightedEdgeSet::const_iterator eit = E.begin (); eit != E.end (); ++ eit )
    f << eit -> first . first << "->" << eit -> first . second << "[color=black,penwidth=5,label=\"" << eit -> second << "\"]" << endl;    

  f << "}" << endl;
  f.close ();

  return;
}

// ConstructDeBruijnGraphSimplified ( "ATTGGTACATTGTGGTACGTACTGACT", 11, 3 );
// ConstructDeBruijnGraphSimplified ( "ATTGATTGACTCTAGATTG", 5, 3 );


