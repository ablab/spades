//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "math/xmath.h"

#include <cstdlib>
#include <set>
#include <unordered_set>

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class VertexPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    VertexPutChecker() { }
    bool Check(VertexId, EdgeId, distance_t) const { return true; }
};

template<class Graph, typename distance_t = size_t>
class EdgeComponentPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const std::set<EdgeId> &edges_;
public:
    EdgeComponentPutChecker(const std::set<EdgeId> &edges) : VertexPutChecker<Graph, distance_t>(), edges_(edges) { }
    bool Check(VertexId, EdgeId edge, distance_t) const {
        return edges_.count(edge) != 0;
    }
};

template<class Graph, typename distance_t = size_t>
class SubgraphPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const std::set<VertexId> &subgraph_;
public:
    SubgraphPutChecker(const std::set<VertexId>& subgraph) : VertexPutChecker<Graph, distance_t>(),
        subgraph_(subgraph) { }
    bool Check(VertexId vertex, EdgeId, distance_t) const {
        return subgraph_.count(vertex) != 0;
    }
};

template<class Graph, typename distance_t = size_t>
class BoundPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const distance_t bound_;

public:
    BoundPutChecker(distance_t bound) :
        bound_(bound) { }
    bool Check(VertexId, EdgeId, distance_t length) const {
        return length <= bound_;
    }
};

template<class Graph, typename distance_t = size_t>
class LengthPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &g_;
    const distance_t length_threshold_;

  public:
    LengthPutChecker(const Graph &g, const distance_t length_threshold) :
        g_(g), length_threshold_(length_threshold) { }
    bool Check(VertexId, EdgeId edge, distance_t ) const {
        return g_.length(edge) <= length_threshold_;
    }
};

template<class Graph, typename distance_t = size_t>
class ForbiddenEdgesPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    std::unordered_set<EdgeId> edges_;
  public:
    ForbiddenEdgesPutChecker(const std::unordered_set<EdgeId> &edges) :
        edges_(edges) { }
    bool Check(VertexId, EdgeId edge, distance_t) const {
        return edges_.find(edge) == edges_.end();
    }
};

template<class Graph, typename distance_t = size_t>
class CoveragePutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    double coverage_low_bound_;
    const Graph &g_;
    const distance_t bound_;
public:
    CoveragePutChecker(double coverage_low_bound, const Graph &g, distance_t bound) : VertexPutChecker<Graph, distance_t>(),
                            coverage_low_bound_(coverage_low_bound), g_(g), bound_(bound){ }
    bool Check(VertexId, EdgeId edge, distance_t length) const {
        return (math::gr(g_.coverage(edge), coverage_low_bound_) && length <= bound_);
    }
};

template<class Tuple, class Graph, typename distance_t = size_t>
class AndPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    Tuple put_checkers_;
  public:
    explicit AndPutChecker(const Tuple &put_checkers) : put_checkers_(put_checkers) { }
    bool Check(VertexId vertex, EdgeId edge, distance_t length) const {
        const auto size = std::tuple_size<Tuple>{};
        return CheckTuple(vertex, edge, length, put_checkers_, std::make_index_sequence<size>{});
    }

  private:
    bool CheckInternal(VertexId, EdgeId, distance_t) const {
        return true;
    }

    template<class CheckerT, class... CheckerTs>
    bool CheckInternal(VertexId vertex, EdgeId edge, distance_t length, const CheckerT &checker,
                       const CheckerTs &... checkers) const {
        return checker.Check(vertex, edge, length) and CheckInternal(vertex, edge, length, checkers...);
    }

    template<size_t... Is>
    bool CheckTuple(VertexId vertex, EdgeId edge, distance_t length, const Tuple &checkers,
                    std::index_sequence<Is...>) const {
        return CheckInternal(vertex, edge, length, std::get<Is>(checkers)...);
    }
};
}
