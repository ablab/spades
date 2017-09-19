//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "func/pred.hpp"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/core/directions.hpp"
#include "assembly_graph/paths/path_finders.hpp"

namespace omnigraph {

template<class Graph>
using EdgePredicate = func::TypedPredicate<typename Graph::EdgeId>;

template<class Graph>
class EdgeCondition : public func::AbstractPredicate<typename Graph::EdgeId> {
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
protected:

    EdgeCondition(const Graph &g)
            : g_(g) {
    }

    const Graph &g() const {
        return g_;
    }

};

template<class Graph>
class IsolatedEdgeCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    bool IsTerminalVertex(VertexId v) const {
        return this->g().IncomingEdgeCount(v) + this->g().OutgoingEdgeCount(v) == 1;
    }

public:
    IsolatedEdgeCondition(const Graph &g) : base(g) {
    }

    bool Check(EdgeId e) const {
        return IsTerminalVertex(this->g().EdgeStart(e)) && IsTerminalVertex(this->g().EdgeEnd(e));
    }

};

template<class Graph>
inline bool HasAlternatives(const Graph &g, typename Graph::EdgeId e) {
    return g.OutgoingEdgeCount(g.EdgeStart(e)) > 1
           && g.IncomingEdgeCount(g.EdgeEnd(e)) > 1;
}


template<class Graph>
class AlternativesPresenceCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

public:

    AlternativesPresenceCondition(const Graph &g)
            : base(g) {

    }

    bool Check(EdgeId e) const {
        return HasAlternatives(this->g(), e);
    }

};

template<class Graph>
func::TypedPredicate<typename Graph::EdgeId> AddAlternativesPresenceCondition(const Graph &g,
                                                                              func::TypedPredicate<typename Graph::EdgeId> condition) {
    return func::And(AlternativesPresenceCondition<Graph>(g), condition);
}


template<class Graph>
class CoverageUpperBound : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef EdgeCondition<Graph> base;
    const double max_coverage_;

public:

    CoverageUpperBound(const Graph &g, double max_coverage)
            : base(g),
              max_coverage_(max_coverage) {
    }

    bool Check(EdgeId e) const {
        return math::le(this->g().coverage(e), max_coverage_);
    }

};

template<class Graph>
class LengthUpperBound : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef EdgeCondition<Graph> base;

    const size_t max_length_;

public:

    LengthUpperBound(const Graph &g, size_t max_length)
            : base(g),
              max_length_(max_length) {
    }

    bool Check(EdgeId e) const {
        return this->g().length(e) <= max_length_;
    }

};

template<class Graph>
class SelfConjugateCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

public:

    SelfConjugateCondition(const Graph& g)
            : base(g) {
    }

    bool Check(EdgeId e) const {
        return e == this->g().conjugate(e);
    }

private:
    DECL_LOGGER("SelfConjugateCondition");
};


}
