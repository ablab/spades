/***************************************************************************
 * Title:          DisjointSet.h 
 * Author:         Mark Chaisson
 * Created:        2007
 * Last modified:  12/18/2008
 *
 * Copyright (c) 2007-2008 The Regents of the University of California
 * All Rights Reserved
 * See file LICENSE for details.
 ***************************************************************************/
#ifndef DISJOINT_SET_H_
#define DISJOINT_SET_H_


// templated disjoint set implementation

namespace DisjointSet {
  template<typename Data_t>
  class DJVertex {
  public:
    DJVertex<Data_t> *parent;
    ssize_t rank;
    Data_t *value;
    DJVertex<Data_t>() {
      parent = this;
      rank = 0;
    }
  };
  
  template<typename Data_t>
  DJVertex<Data_t> *Find(DJVertex<Data_t> *node) {
    if (node != node->parent )
      node->parent = Find(node->parent);
    return node->parent;
  }

  template<typename Data_t>
  void Union(DJVertex<Data_t> *x, DJVertex<Data_t> *y) {
    Link(Find(x), Find(y));
  }

  template<typename Data_t>
  void Link(DJVertex<Data_t> *x, DJVertex<Data_t> *y) {
    if (x->rank > y->rank) 
      y->parent = x;
    else {
      x->parent = y;
      if (x->rank == y->rank)
				y->rank++;
    }
  }


  template<typename Vertex_t>
  void CreateForest(std::vector<Vertex_t> &vertices, std::vector<DJVertex<Vertex_t> > &forest) {
    forest.resize(vertices.size());
    ssize_t i;
    for (i = 0; i < vertices.size(); i++) {
      forest[i].value = &(vertices[i]);
      forest[i].parent = &(forest[i]);
    }
  }
};


#endif
