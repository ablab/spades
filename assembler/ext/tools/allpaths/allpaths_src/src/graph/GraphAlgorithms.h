///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef GRAPH__GRAPH_ALGORITHMS__H_
#define GRAPH__GRAPH_ALGORITHMS__H_



class UnionFind
{
  vec<size_t> iv_parent;   // index of parent vertex for each vertex
  vec<size_t> rank;       // to keep tree mergers balanced

public:
  UnionFind(const size_t nv) : iv_parent(nv), rank(nv, 0) 
  {
    for (size_t iv = 0; iv < nv; iv++)
      iv_parent[iv] = iv;  // in the begining every node is its own tree in the forest
  }

  size_t n_vertexes() const { return iv_parent.size(); }

  size_t iv_root(const size_t iv) const // recursively find root of tree
  {
    const size_t ivp = iv_parent[iv];
    return (ivp == iv) ? ivp : iv_root(ivp);
  }

  size_t iv_root(const size_t iv)
  {
    const size_t ivp = iv_parent[iv];
    return (ivp == iv) ? 
      ivp : 
      iv_parent[iv] = iv_root(ivp); // path-compression (non-const method)
  }
  
  void unite(const size_t iv0, const size_t iv1)
  {
    const size_t iv0_root = iv_root(iv0);
    const size_t iv1_root = iv_root(iv1);
    
    // union-by-rank to keep tree balanced

    if (iv0_root != iv1_root) {
      if      (rank[iv0_root] < rank[iv1_root])  iv_parent[iv0_root] = iv1_root;
      else if (rank[iv0_root] > rank[iv1_root])  iv_parent[iv1_root] = iv0_root;
      else {
        iv_parent[iv1_root] = iv0_root;
        rank[iv0_root]++;
      }
    }
  }

  void report_print()
  {
    map<size_t, size_t> nv_tree;
    const size_t nv = n_vertexes();
    for (size_t iv = 0; iv < nv; iv++)
      nv_tree[iv_root(iv)]++;

    cout << "n_trees= " << nv_tree.size() << endl;

    for (map<size_t, size_t>::const_iterator it = nv_tree.begin();
         it != nv_tree.end(); it++)
      cout << "n_vertexes[" << setw(6) << (*it).first << "]= " << setw(6) << (*it).second << endl;
  }
  
};










template<class WEIGHT_t>
class EdgeWeight
{
public:
  size_t iv0;        // vertex index 0
  size_t iv1;        // vertex index 1
  WEIGHT_t weight;   // edge weight
  
  EdgeWeight(const size_t _iv0, const size_t _iv1, const WEIGHT_t & _w)
    : iv0(_iv0), iv1(_iv1), weight(_w) {}
};









template<class WEIGHT_t>
void minimum_spanning_tree_kruskal(const vec<EdgeWeight<WEIGHT_t> > & edges,
                                   UnionFind                        * tree_p,
                                   vec<size_t>                      * i_edges_p)
{
  // ---- sort edges by weight
  
  vec<pair<WEIGHT_t, size_t> > weight_ie;
  
  const size_t ne = edges.size();
  
  for (size_t ie = 0; ie < ne; ie++) 
    weight_ie.push_back(make_pair(edges[ie].weight, ie));
  
  sort(weight_ie.begin(), weight_ie.end());
  

  // ---- start picking edges from the lightest one

  for (size_t iie = 0; iie < ne; iie++) {
    const size_t ie = weight_ie[iie].second;
    const EdgeWeight<WEIGHT_t> & edge = edges[ie];
    if (tree_p->iv_root(edge.iv0) != tree_p->iv_root(edge.iv1)) {
      tree_p->unite(edge.iv0, edge.iv1);
      i_edges_p->push_back(ie);
    }
  }
}
                                  

template<class WEIGHT_t>
void minimum_spanning_tree_kruskal(const size_t                       n_vertexes,
                                   const vec<EdgeWeight<WEIGHT_t> > & edges,
                                   vec<size_t>                      * i_edges_p)
{
  UnionFind tree(n_vertexes);
  minimum_spanning_tree_kruskal(edges, & tree, i_edges_p);
}


















#endif
