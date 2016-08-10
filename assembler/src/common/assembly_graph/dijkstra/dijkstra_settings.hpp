//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "length_calculator.hpp"
#include "vertex_process_checker.hpp"
#include "vertex_put_checker.hpp"
#include "neighbours_iterator.hpp"

namespace omnigraph {

template<class Graph,
        class LengthCalculator,
        class VertexProcessChecker,
        class VertexPutChecker,
        class NeighbourIteratorFactory,
        typename distance_t = size_t>
class ComposedDijkstraSettings {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    LengthCalculator len_calc_;
    VertexProcessChecker vert_proc_checker_;
    VertexPutChecker vert_put_checker_;
    NeighbourIteratorFactory neigh_iter_factory_;

public:
    typedef LengthCalculator LC;
    typedef VertexProcessChecker VPrC;
    typedef VertexPutChecker VPuC;
    typedef NeighbourIteratorFactory NIF;

    ComposedDijkstraSettings(LengthCalculator len_calc,
            VertexProcessChecker vert_proc_checker,
            VertexPutChecker vert_put_checker,
            NeighbourIteratorFactory neigh_iter_factory) :
                len_calc_(len_calc),
                vert_proc_checker_(vert_proc_checker),
                vert_put_checker_(vert_put_checker),
                neigh_iter_factory_(neigh_iter_factory) { }

    void Init(VertexId /*vertex*/){
    }

    distance_t GetLength(EdgeId edge) const{
        return len_calc_.GetLength(edge);
    }

    bool CheckProcessVertex(VertexId vertex, distance_t distance){
        return vert_proc_checker_.Check(vertex, distance);
    }

    bool CheckPutVertex(VertexId vertex, EdgeId edge, distance_t length) const{
        return vert_put_checker_.Check(vertex, edge, length);
    }

    typename NeighbourIteratorFactory::NeighbourIterator GetIterator(VertexId vertex) {
        return neigh_iter_factory_.CreateIterator(vertex);
    }
};

template<class Graph, class NeighbourIteratorFactory, typename distance_t = size_t>
class CountingDijkstraSettings {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;

    NeighbourIteratorFactory neigh_iter_factory_;
    static const distance_t inf = 100000000;
    const size_t max_size_;
    const size_t edge_length_bound_;
    mutable size_t current_;

public:
    CountingDijkstraSettings(const Graph &graph,
            NeighbourIteratorFactory neigh_iter_factory,
            size_t max_size, size_t edge_length_bound) :
           graph_(graph),
           neigh_iter_factory_(neigh_iter_factory),
        max_size_(max_size),
        edge_length_bound_(edge_length_bound),
        current_(0) { }

    void Init(VertexId /*vertex*/){
        current_ = 0;
    }

    distance_t GetLength(EdgeId edge) const{
        if (graph_.length(edge) <= edge_length_bound_)
            return graph_.length(edge);
        return inf;
    }

    bool CheckProcessVertex(VertexId , distance_t ){
        return current_ < max_size_;
    }

    bool CheckPutVertex(VertexId , EdgeId edge, distance_t ) const{
        if (current_ < max_size_)
            ++current_;
        if (current_ < max_size_ && GetLength(edge) < inf)
            return true;
        return false;
    }

    typename NeighbourIteratorFactory::NeighbourIterator GetIterator(VertexId vertex) {
        return neigh_iter_factory_.CreateIterator(vertex);
    }
};

}
