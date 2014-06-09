//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "func.hpp"
#include "omni_utils.hpp"

namespace omnigraph {

using namespace func;

template<class Graph>
class EdgeCondition : public Predicate<typename Graph::EdgeId> {
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
 protected:

    EdgeCondition(const Graph& g)
            : g_(g) {
    }

    const Graph& g() const {
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
	IsolatedEdgeCondition(const Graph& g) : base(g) {
	}

    bool Check(EdgeId e) const {
        return IsTerminalVertex(this->g().EdgeStart(e)) && IsTerminalVertex(this->g().EdgeEnd(e));
    }

};

template<class Graph>
class AlternativesPresenceCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

 public:

    AlternativesPresenceCondition(const Graph& g)
            : base(g) {

    }

    bool Check(EdgeId e) const {
        return this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) > 1
                && this->g().IncomingEdgeCount(this->g().EdgeEnd(e)) > 1;
    }

};

template<class Graph>
class CoverageUpperBound : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef EdgeCondition<Graph> base;
    const double max_coverage_;

 public:

    CoverageUpperBound(const Graph& g, double max_coverage)
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

    LengthUpperBound(const Graph& g, size_t max_length)
            : base(g),
              max_length_(max_length) {
    }

    bool Check(EdgeId e) const {
        return this->g().length(e) <= max_length_;
    }

};

template<class Graph, class PathFinder>
class PathLengthLowerBound : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    PathFinder path_finder_;
    size_t min_length_;

    ForwardDirection<Graph> forward_;
    BackwardDirection<Graph> backward_;

    size_t CumulativePathLength(EdgeId e, const AbstractDirection<Graph>& direction) const {
        return CumulativeLength(this->g(), path_finder_(e, direction));
    }

 public:
    PathLengthLowerBound(const Graph& g, const PathFinder& path_finder,
                         size_t min_length)
            : base(g),
              path_finder_(path_finder),
              min_length_(min_length),
              forward_(g),
              backward_(g) {

    }

    bool Check(EdgeId e) const {
        size_t forward = CumulativePathLength(e, forward_);
        size_t backward = CumulativePathLength(e, backward_);
        //checking that path was trivial in one of directions
        VERIFY(forward == this->g().length(e) || backward == this->g().length(e));
        return std::max(forward, backward) >= min_length_;
    }
};

template<class Graph, class PathFinder>
std::shared_ptr<PathLengthLowerBound<Graph, PathFinder> >
MakePathLengthLowerBound(const Graph& g, const PathFinder& path_finder, size_t min_length) {
    return std::make_shared<PathLengthLowerBound<Graph, PathFinder>>(g, path_finder,
                                                                min_length);
}

template<class Graph>
class UniquenessPlausabilityCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    virtual bool CheckUniqueness(EdgeId e, bool forward) const = 0;

    virtual bool CheckPlausibility(EdgeId e, bool forward) const = 0;

    bool SingleUnique(const vector<EdgeId>& edges, bool forward) const {
        return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
    }

    bool ExistPlausible(EdgeId init_e, const vector<EdgeId>& edges,
                        bool forward) const {
        FOREACH(EdgeId e, edges) {
            if (e == init_e)
                continue;
            if (CheckPlausibility(e, forward)) {
                return true;
            }
        }
        return false;
    }

    bool Check(EdgeId e, const AbstractDirection<Graph>& direction) const {
        return SingleUnique(direction.IncomingEdges(direction.EdgeStart(e)),
                            !direction.IsForward())
                && ExistPlausible(
                        e, direction.OutgoingEdges(direction.EdgeStart(e)),
                        direction.IsForward());
    }

 public:

    UniquenessPlausabilityCondition(const Graph& g)
            : base(g) {

    }

    bool Check(EdgeId e) const {
        return Check(e, ForwardDirection<Graph>(this->g()))
                || Check(e, BackwardDirection<Graph>(this->g()));
    }

};

template<class Graph>
class PredicateUniquenessPlausabilityCondition :
        public UniquenessPlausabilityCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef shared_ptr<Predicate<EdgeId>> EdgePredicate;
    typedef UniquenessPlausabilityCondition<Graph> base;

    EdgePredicate uniqueness_condition_;
    EdgePredicate plausiblity_condition_;

    bool CheckUniqueness(EdgeId e, bool) const {
        return uniqueness_condition_->Check(e);
    }

    bool CheckPlausibility(EdgeId e, bool) const {
        return plausiblity_condition_->Check(e);
    }

 public:

    PredicateUniquenessPlausabilityCondition(
            const Graph& g, EdgePredicate uniqueness_condition,
            EdgePredicate plausiblity_condition)
            : base(g),
              uniqueness_condition_(uniqueness_condition),
              plausiblity_condition_(plausiblity_condition) {
    }

};

template<class Graph>
class DefaultUniquenessPlausabilityCondition :
        public PredicateUniquenessPlausabilityCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef shared_ptr<Predicate<EdgeId>> EdgePredicate;
    typedef PredicateUniquenessPlausabilityCondition<Graph> base;

 public:

    DefaultUniquenessPlausabilityCondition(const Graph& g,
                                           size_t uniqueness_length,
                                           size_t plausibility_length)
            : base(g,
                   MakePathLengthLowerBound(g, UniquePathFinder<Graph>(g),
                                            uniqueness_length),
                   MakePathLengthLowerBound(
                           g,
                           PlausiblePathFinder<Graph>(g,
                                                      2 * plausibility_length),
                           plausibility_length)) {
    }

};

}
