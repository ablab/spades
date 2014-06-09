//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "standard_base.hpp"
#include "simple_tools.hpp"
#include "xmath.h"

#include "omni/action_handlers.hpp"
#include "omni/graph_iterators.hpp"
#include "omni/mapping_path.hpp"

#include <cmath>
#include <ctime>

namespace omnigraph {

template<class Graph>
struct CoverageComparator {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& graph_;
 public:
    CoverageComparator(const Graph &graph)
            : graph_(graph) {
    }

    /**
     * Standard comparator function as used in collections.
     */
    bool operator()(EdgeId edge1, EdgeId edge2) const {
        if (math::eq(graph_.coverage(edge1), graph_.coverage(edge2))) {
            return edge1 < edge2;
        }
        return math::ls(graph_.coverage(edge1), graph_.coverage(edge2));
    }
};

/**
 * This class defines which edge is more likely to be tip. In this case we just assume shorter edges
 * are more likely tips then longer ones.
 */
template<class Graph>
struct LengthComparator {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph& graph_;
 public:
    /**
     * TipComparator should never be created with default constructor but it is necessary on order for
     * code to compile.
     */
    //  TipComparator() {
    //    VERIFY(false);
    //  }
    /**
     * Construct TipComparator for given graph
     * @param graph graph for which comparator is created
     */
    LengthComparator(const Graph &graph)
            : graph_(graph) {
    }

    /**
     * Standard comparator function as used in collections.
     */
    bool operator()(EdgeId edge1, EdgeId edge2) const {
        if (graph_.length(edge1) == graph_.length(edge2)) {
            return edge1 < edge2;
        }
        return graph_.length(edge1) < graph_.length(edge2);
    }
};

template<class Graph>
size_t CumulativeLength(const Graph& g,
                        const std::vector<typename Graph::EdgeId>& path) {
    size_t s = 0;
    for (auto it = path.begin(); it != path.end(); ++it)
        s += g.length(*it);

    return s;
}

template<class Graph>
double AvgCoverage(const Graph& g,
                   const std::vector<typename Graph::EdgeId>& path) {
    double unnormalized_coverage = 0;
    size_t path_length = 0;
    for (auto edge : path) {
        size_t length = g.length(edge);
        path_length += length;
        unnormalized_coverage += g.coverage(edge) * (double) length;
    }
    return unnormalized_coverage / (double) path_length;
}

template<class Graph>
class AbstractDirection {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;

 protected:
    const Graph &graph() const {
        return graph_;
    }

 public:
    AbstractDirection(const Graph& graph)
            : graph_(graph) {}

    virtual ~AbstractDirection() {}

    virtual const std::vector<EdgeId> OutgoingEdges(VertexId v) const = 0;
    virtual const std::vector<EdgeId> IncomingEdges(VertexId v) const = 0;

    virtual size_t OutgoingEdgeCount(VertexId v) const = 0;
    virtual size_t IncomingEdgeCount(VertexId v) const = 0;

    virtual VertexId EdgeStart(EdgeId edge) const = 0;
    virtual VertexId EdgeEnd(EdgeId edge) const = 0;

    bool CheckUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueOutgoingEdge(VertexId v) const {
        return OutgoingEdges(v)[0];
    }

    bool CheckUniqueIncomingEdge(VertexId v) const {
        return IncomingEdgeCount(v) == 1;
    }

    EdgeId GetUniqueIncomingEdge(VertexId v) const {
        return IncomingEdges(v)[0];
    }

    virtual bool IsForward() const = 0;
};

template<class Graph>
class ForwardDirection : public AbstractDirection<Graph> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
 public:
    ForwardDirection(const Graph &graph)
            : AbstractDirection<Graph>(graph) {
    }

    virtual const std::vector<EdgeId> OutgoingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().out_begin(v), this->graph().out_end(v));
    }

    virtual const std::vector<EdgeId> IncomingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().in_begin(v), this->graph().in_end(v));
    }

    virtual size_t OutgoingEdgeCount(VertexId v) const {
        return this->graph().OutgoingEdgeCount(v);
    }

    virtual size_t IncomingEdgeCount(VertexId v) const {
        return this->graph().IncomingEdgeCount(v);
    }

    virtual VertexId EdgeStart(EdgeId edge) const {
        return this->graph().EdgeStart(edge);
    }

    virtual VertexId EdgeEnd(EdgeId edge) const {
        return this->graph().EdgeEnd(edge);
    }

    bool IsForward() const {
        return true;
    }
};

template<class Graph>
class BackwardDirection : public AbstractDirection<Graph> {
 private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
 public:
    BackwardDirection(const Graph &graph)
            : AbstractDirection<Graph>(graph) {
    }

    virtual const std::vector<EdgeId> OutgoingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().in_begin(v), this->graph().in_end(v));
    }

    virtual const std::vector<EdgeId> IncomingEdges(VertexId v) const {
        return std::vector<EdgeId>(this->graph().out_begin(v), this->graph().out_end(v));
    }

    virtual size_t OutgoingEdgeCount(VertexId v) const {
        return this->graph().IncomingEdgeCount(v);
    }

    virtual size_t IncomingEdgeCount(VertexId v) const {
        return this->graph().OutgoingEdgeCount(v);
    }

    virtual VertexId EdgeStart(EdgeId edge) const {
        return this->graph().EdgeEnd(edge);
    }

    virtual VertexId EdgeEnd(EdgeId edge) const {
        return this->graph().EdgeStart(edge);
    }

    bool IsForward() const {
        return false;
    }

};

template<class Graph>
class UniquePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;
 public:
    //todo use length bound if needed
    UniquePathFinder(const Graph& graph, size_t /*length_bound*/ =
                             std::numeric_limits<size_t>::max())
            : graph_(graph) {}

    std::vector<EdgeId> operator()(EdgeId e,
                                   const AbstractDirection<Graph> &direction) const {
        std::vector<EdgeId> answer;
        EdgeId curr = e;
        answer.push_back(curr);
        std::set<EdgeId> was;
        while (direction.CheckUniqueOutgoingEdge(direction.EdgeEnd(curr))) {
            curr = direction.GetUniqueOutgoingEdge(direction.EdgeEnd(curr));
            if (was.count(curr) > 0)
                break;
            was.insert(curr);
            answer.push_back(curr);
        }
        return answer;
    }

    std::vector<EdgeId> UniquePathForward(EdgeId e) const {
        return this->operator()(e, ForwardDirection<Graph>(graph_));
    }

    std::vector<EdgeId> UniquePathBackward(EdgeId e) const {
        auto tmp = this->operator()(e, BackwardDirection<Graph>(graph_));
        return std::vector<EdgeId>(tmp.rbegin(), tmp.rend());
    }

};

template<class Graph>
class TrivialPathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

 public:
    TrivialPathFinder(const Graph&, size_t = 0) {}

    std::vector<EdgeId> operator()(EdgeId e, const AbstractDirection<Graph> &) const {
        return {e};
    }

};

template<class Graph>
class PlausiblePathFinder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    //todo remove graph_ field???
    const Graph& graph_;
    const size_t length_bound_;

    class DFS {
     private:
        const Graph &graph_;
        const AbstractDirection<Graph> &direction_;
        const size_t length_bound_;

        std::pair<size_t, EdgeId> find(EdgeId edge, size_t length) {
            length += graph_.length(edge);
            VertexId cross = direction_.EdgeEnd(edge);
            auto result = make_pair(length, edge);
            if (length < length_bound_
                    && direction_.CheckUniqueIncomingEdge(cross)) {
                std::vector<EdgeId> outgoing = direction_.OutgoingEdges(cross);
                for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
                    auto candidate = find(*it, length);
                    if (candidate.first > result.first)
                        result = candidate;
                }
            }
            return result;
        }

        std::vector<EdgeId> RestoreAnswer(EdgeId start, EdgeId end) {
            std::vector<EdgeId> result;
            while (end != start) {
                result.push_back(end);
                end = direction_.GetUniqueIncomingEdge(direction_.EdgeStart(end));
            }
            result.push_back(start);
            return std::vector<EdgeId>(result.rbegin(), result.rend());
        }

     public:
        DFS(const Graph &graph, const AbstractDirection<Graph> &direction,
            size_t length_bound)
                : graph_(graph),
                  direction_(direction),
                  length_bound_(length_bound) {
        }

        std::vector<EdgeId> find(EdgeId edge) {
            return RestoreAnswer(edge, find(edge, 0).second);
        }
    };

 public:
    PlausiblePathFinder(const Graph& graph, size_t length_bound)
            : graph_(graph),
              length_bound_(length_bound) {}

    std::vector<EdgeId> operator()(EdgeId e,
                                   const AbstractDirection<Graph> &direction) const {
        return DFS(graph_, direction, length_bound_).find(e);
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
        FOREACH (EdgeId in_e, graph_.IncomingEdges(a)) {
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
class DominatedSetFinder {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    VertexId start_vertex_;
    size_t max_length_;
    size_t max_count_;

    size_t cnt_;
    std::map<VertexId, Range> dominated_;

    bool CheckCanBeProcessed(VertexId v) const {
        DEBUG( "Check if vertex " << g_.str(v) << " is dominated close neighbour");
        FOREACH (EdgeId e, g_.IncomingEdges(v)) {
            if (dominated_.count(g_.EdgeStart(e)) == 0) {
                DEBUG( "Blocked by external vertex " << g_.int_id(g_.EdgeStart(e)) << " that starts edge " << g_.int_id(e));
                DEBUG("Check fail");
                return false;
            }
        }
        DEBUG("Check ok");
        return true;
    }

    void UpdateCanBeProcessed(VertexId v,
                              std::queue<VertexId>& can_be_processed) const {
        DEBUG("Updating can be processed")
        FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
            DEBUG("Considering edge " << ToString(e));
            VertexId neighbour_v = g_.EdgeEnd(e);
            if (CheckCanBeProcessed(neighbour_v)) {
                can_be_processed.push(neighbour_v);
            }
        }
    }

    Range NeighbourDistanceRange(VertexId v, bool dominated_only = true) const {
        DEBUG("Counting distance range for vertex " << g_.str(v));
        size_t min = numeric_limits<size_t>::max();
        size_t max = 0;
        VERIFY(g_.IncomingEdgeCount(v) > 0);
        VERIFY(!dominated_only || CheckCanBeProcessed(v));
        FOREACH (EdgeId e, g_.IncomingEdges(v)) {
            //in case of dominated_only == false
            if (dominated_.count(g_.EdgeStart(e)) == 0)
                continue;
            Range range = dominated_.find(g_.EdgeStart(e))->second;
            range.shift((int) g_.length(e));
            DEBUG("Edge " << g_.str(e) << " provide distance range " << range);
            if (range.start_pos < min)
                min = range.start_pos;
            if (range.end_pos > max)
                max = range.end_pos;
        }
        VERIFY((max > 0) && (min < numeric_limits<size_t>::max()) && (min <= max));
        Range answer(min, max);
        DEBUG("Range " << answer);
        return answer;
    }

    bool CheckNoEdgeToStart(VertexId v) {
        FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
            if (g_.EdgeEnd(e) == start_vertex_) {
                return false;
            }
        }
        return true;
    }

 public:
    DominatedSetFinder(const Graph& g, VertexId v, size_t max_length = -1ul,
                       size_t max_count = -1ul)
            : g_(g),
              start_vertex_(v),
              max_length_(max_length),
              max_count_(max_count),
              cnt_(0) {

    }

    //true if no thresholds exceeded
    bool FillDominated() {
        DEBUG("Adding starting vertex " << g_.str(start_vertex_) << " to dominated set");
        dominated_.insert(make_pair(start_vertex_, Range(0, 0)));
        cnt_++;
        std::queue<VertexId> can_be_processed;
        UpdateCanBeProcessed(start_vertex_, can_be_processed);
        while (!can_be_processed.empty()) {
            if (++cnt_ > max_count_) {
                return false;
            }
            VertexId v = can_be_processed.front();
            can_be_processed.pop();
            Range r = NeighbourDistanceRange(v);
            if (r.start_pos > max_length_) {
                return false;
            }
            //Currently dominated vertices cannot have edge to start vertex
            if (CheckNoEdgeToStart(v)) {
                DEBUG("Adding vertex " << g_.str(v) << " to dominated set");
                dominated_.insert(make_pair(v, r));
                UpdateCanBeProcessed(v, can_be_processed);
            }
        }
        return true;
    }

    const map<VertexId, Range>& dominated() const {
        return dominated_;
    }

    //little meaning if FillDominated returned false
    const map<VertexId, Range> CountBorder() const {
        map<VertexId, Range> border;
        FOREACH(VertexId v, key_set(border)) {
            FOREACH(EdgeId e, g_.OutgoingEdges(v)) {
                VertexId e_end = g_.EdgeEnd(e);
                if (dominated_.count(e_end) == 0) {
                    border[e_end] = NeighbourDistanceRange(e_end, false);
                }
            }
        }
        return border;
    }

};

inline size_t PairInfoPathLengthUpperBound(size_t k, size_t insert_size,
                                           double delta) {
    double answer = 0. + (double) insert_size + delta - (double) k - 2.;
    VERIFY(math::gr(answer, 0.));
    return (size_t)std::floor(answer);
}

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l1, size_t l2,
                                           int gap, double delta) {
    double answer = 0. + (double) gap + (double) k + 2. - (double) l1 - (double) l2 - delta;
    return math::gr(answer, 0.) ? (size_t)std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */
