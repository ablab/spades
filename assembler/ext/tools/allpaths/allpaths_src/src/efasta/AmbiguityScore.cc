///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file AmbiguityScore.cc
 * \author tsharpe
 * \date Dec 14, 2011
 *
 * \brief
 */
#include "efasta/AmbiguityScore.h"
#include "Alignment.h"
#include "pairwise_aligners/SmithWatFree.h"
#include <utility>

#include "graph/GraphAlgorithms.h"




int AmbiguityScore(const vec<basevector>& x)
{
  const size_t nv = x.size();
  vec<EdgeWeight<int> > edges;
  
  {
    int best_loc;
    alignment a;
    for (size_t iv = 0; iv < nv; iv++) { 
      for (size_t jv = iv + 1; jv < nv; jv++) {
        
        int score;
        if (x[iv].size() == 0 || x[jv].size() == 0)
          score = Abs(x[iv].isize() - x[jv].isize());
        else if (x[iv].size() <= x[jv].size())
          score = SmithWatFree(x[iv], x[jv], best_loc, a, true, true, 1, 1, 1);
        else
          score = SmithWatFree(x[jv], x[iv], best_loc, a, true, true, 1, 1, 1);

        edges.push_back(EdgeWeight<int>(iv, jv, score));
      }
    }
  }


  vec<size_t> i_edges_tree;
  minimum_spanning_tree_kruskal(nv, edges, & i_edges_tree);
  
  const size_t ne = i_edges_tree.size();

  int total = nv - 2;
  for (size_t iie = 0; iie < ne; iie++) {
    const size_t ie = i_edges_tree[iie];
    total += edges[ie].weight;
  }

  return total; 
}






int64_t AmbiguityScoreCost( const vec<basevector>& x )
{    int64_t cost = 0;
     for ( int i = 0; i < x.isize( ); i++ )
     {    for ( int j = i+1; j < x.isize( ); j++ )
               cost += x[i].size( ) * x[j].size( );    }
     return cost;    }
