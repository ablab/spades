#include "toyexamples.h"

#include <string>
#include <set>
#include <map>
#include <iostream>
#include <cassert>

#include "../graphVisualizer.hpp"

using namespace std;
using namespace gvis;

typedef map < string, unsigned > CKMerSet;

extern void ConstructDeBruijnGraph ( string genome, unsigned k )//, string filename )
{
  //ofstream & ofs

  unsigned i = 0;// just counter
  unsigned total_kmer_num = 0;

  GraphScheme < unsigned > G ( "DEBRUIJNGRAPH" );

  CKMerSet S;

  // adding vertices
  for ( i = 0; i != genome.size () - k + 1; ++i )
    if ( S.find ( genome.substr ( i, k ) ) == S.end () )
    {
      ++total_kmer_num;
      S.insert ( make_pair ( genome.substr ( i, k ), total_kmer_num ) );
      G.addVertex ( total_kmer_num, genome.substr ( i, k ) );
      //debug
      //cout << "\n adding vertex " << total_kmer_num << ": " << genome.substr ( i, k ) << flush;
    }

  // adding edges
  for ( i = 0; i != genome.size () - ( k + 1 )+ 1; ++i )
  {
    //string edgelabel = genome.substr ( i, k + 1);
    CKMerSet::const_iterator from = S.find ( genome.substr ( i,     k ) );
    CKMerSet::const_iterator to   = S.find ( genome.substr ( i + 1, k ) );
    assert ( from != S.end () && to != S.end () );
    G.addEdge ( from -> second, to -> second, genome.substr ( i, k + 1 ) );
  }

  G.output ();

  return;
}
