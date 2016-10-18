#pragma once

#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "assembly_graph/core/directions.hpp"

namespace omnigraph {

template<class Graph, class PathFinder>
class PathLengthLowerBound : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    PathFinder path_finder_;
    size_t min_length_;

    ForwardDirection<Graph> forward_;
    BackwardDirection<Graph> backward_;

    size_t CumulativePathLength(EdgeId e, const AbstractDirection<Graph> &direction) const {
        return CumulativeLength(this->g(), path_finder_(e, direction));
    }

public:
    PathLengthLowerBound(const Graph &g, const PathFinder &path_finder,
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
EdgePredicate<Graph>
MakePathLengthLowerBound(const Graph &g, const PathFinder &path_finder, size_t min_length) {
    return PathLengthLowerBound<Graph, PathFinder>(g, path_finder, min_length);
}

template<class Graph>
EdgePredicate<Graph>
UniquePathLengthLowerBound(const Graph &g, size_t min_length) {
    return MakePathLengthLowerBound(g, UniquePathFinder<Graph>(g), min_length);
}

template<class Graph>
EdgePredicate<Graph>
UniqueIncomingPathLengthLowerBound(const Graph &g, size_t min_length) {
    return [&] (typename Graph::EdgeId e) {
        typename Graph::VertexId v = g.EdgeStart(e);
        return g.CheckUniqueIncomingEdge(v) &&
                UniquePathLengthLowerBound(g, min_length)(g.GetUniqueIncomingEdge(v));
    };
}

//todo can disconnect uniqueness and plausibility conditions, since graph is always conjugate!
template<class Graph>
class UniquenessPlausabilityCondition : public EdgeCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef EdgeCondition<Graph> base;

    virtual bool CheckUniqueness(EdgeId e, bool forward) const = 0;

    virtual bool CheckPlausibility(EdgeId e, bool forward) const = 0;

    bool SingleUnique(const vector<EdgeId> &edges, bool forward) const {
        return edges.size() == 1 && CheckUniqueness(*edges.begin(), forward);
    }

    bool ExistPlausible(EdgeId init_e, const vector<EdgeId> &edges,
                        bool forward) const {
        for (EdgeId e : edges) {
            if (e == init_e)
                continue;
            if (CheckPlausibility(e, forward)) {
                return true;
            }
        }
        return false;
    }

    bool Check(EdgeId e, const AbstractDirection<Graph> &direction) const {
        return SingleUnique(direction.IncomingEdges(direction.EdgeStart(e)),
                            !direction.IsForward())
               && ExistPlausible(
                e, direction.OutgoingEdges(direction.EdgeStart(e)),
                direction.IsForward());
    }

public:

    UniquenessPlausabilityCondition(const Graph &g)
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
    typedef UniquenessPlausabilityCondition<Graph> base;

    EdgePredicate<Graph> uniqueness_condition_;
    EdgePredicate<Graph> plausiblity_condition_;

    bool CheckUniqueness(EdgeId e, bool) const {
        return uniqueness_condition_(e);
    }

    bool CheckPlausibility(EdgeId e, bool) const {
        return plausiblity_condition_(e);
    }

public:

    PredicateUniquenessPlausabilityCondition(
            const Graph &g, EdgePredicate<Graph> uniqueness_condition,
            EdgePredicate<Graph> plausiblity_condition)
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
    typedef PredicateUniquenessPlausabilityCondition<Graph> base;

public:

    DefaultUniquenessPlausabilityCondition(const Graph &g,
                                           size_t uniqueness_length,
                                           size_t plausibility_length)
            : base(g,
                   UniquePathLengthLowerBound(g, uniqueness_length),
                   MakePathLengthLowerBound(g,
                                            PlausiblePathFinder<Graph>(g, 2 * plausibility_length),
                                            plausibility_length)) {
    }

};

template<class Graph>
class MultiplicityCounter {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
    size_t uniqueness_length_;
    size_t max_depth_;

    bool search(VertexId a, VertexId start, EdgeId e, size_t depth,
                std::set<VertexId> &was, pair<size_t, size_t> &result) const {
        if (depth > max_depth_)
            return false;
        if (was.count(a) == 1)
            return true;
        was.insert(a);
        if (graph_.OutgoingEdgeCount(a) == 0
            || graph_.IncomingEdgeCount(a) == 0)
            return false;
        for (auto I = graph_.out_begin(a), E = graph_.out_end(a); I != E; ++I) {
            if (*I == e) {
                if (a != start) {
                    return false;
                }
            } else {
                if (graph_.length(*I) >= uniqueness_length_) {
                    result.second++;
                } else {
                    if (!search(graph_.EdgeEnd(*I), start, e,
                                depth + 1 /*graph_.length(*it)*/, was, result))
                        return false;
                }
            }
        }
        for (EdgeId in_e : graph_.IncomingEdges(a)) {
            if (in_e == e) {
                if (a != start) {
                    return false;
                }
            } else {
                if (graph_.length(in_e) >= uniqueness_length_) {
                    result.first++;
                } else {
                    if (!search(graph_.EdgeStart(in_e), start, e,
                                depth + 1 /*graph_.length(*it)*/, was, result))
                        return false;
                }
            }
        }
        return true;
    }

public:
    MultiplicityCounter(const Graph &graph, size_t uniqueness_length,
                        size_t max_depth)
            : graph_(graph),
              uniqueness_length_(uniqueness_length),
              max_depth_(max_depth) {
    }

    size_t count(EdgeId e, VertexId start) const {
        std::pair<size_t, size_t> result;
        std::set<VertexId> was;
        bool valid = search(start, start, e, 0, was, result);
        if (!valid) {
            return (size_t) (-1);
        }
        if (graph_.EdgeStart(e) == start) {
            if (result.first < result.second) {
                return (size_t) (-1);
            }
            return result.first - result.second;
        } else {
            if (result.first > result.second) {
                return (size_t) (-1);
            }
            return -result.first + result.second;
        }
    }
};

template<class Graph>
class MultiplicityCountingCondition : public UniquenessPlausabilityCondition<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef UniquenessPlausabilityCondition<Graph> base;

    MultiplicityCounter<Graph> multiplicity_counter_;
    EdgePredicate<Graph> plausiblity_condition_;

public:
    bool CheckUniqueness(EdgeId e, bool forward) const {
        TRACE( "Checking " << this->g().int_id(e) << " for uniqueness in " << (forward ? "forward" : "backward") << " direction");
        VertexId start =
                forward ? this->g().EdgeEnd(e) : this->g().EdgeStart(e);
        bool result = multiplicity_counter_.count(e, start) <= 1;
        TRACE( "Edge " << this->g().int_id(e) << " is" << (result ? "" : " not") << " unique");
        return result;
    }

    bool CheckPlausibility(EdgeId e, bool) const {
        return plausiblity_condition_(e);
    }

    MultiplicityCountingCondition(const Graph& g, size_t uniqueness_length,
                                  EdgePredicate<Graph> plausiblity_condition)
            :
    //todo why 8???
            base(g),
            multiplicity_counter_(g, uniqueness_length, 8),
            plausiblity_condition_(plausiblity_condition) {

    }

private:

    DECL_LOGGER("MultiplicityCountingCondition");
};


}
