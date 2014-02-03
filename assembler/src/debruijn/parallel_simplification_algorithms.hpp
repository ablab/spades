#pragma once

#include "standard_base.hpp"
#include "omni/graph_processing_algorithm.hpp"
#include "omni/basic_edge_conditions.hpp"
#include "omni/bulge_remover.hpp"
#include <boost/range/adaptor/reversed.hpp>

namespace debruijn {

namespace simplification {

template<class Graph>
class ParallelTipClippingFunctor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef boost::function<void(EdgeId)> HandlerF;

    Graph& g_;
    size_t length_bound_;
    HandlerF handler_f_;

    size_t LockingIncomingCount(VertexId v) const {
        size_t answer;
        v->Lock();
        answer = g_.IncomingEdgeCount(v);
        v->Unlock();
        return answer;
    }

    size_t LockingOutgoingCount(VertexId v) const {
        size_t answer;
        v->Lock();
        answer = g_.OutgoingEdgeCount(v);
        v->Unlock();
        return answer;
    }

    bool IsIncomingTip(EdgeId e) const {
        return g_.length(e) <= length_bound_ && LockingIncomingEdgeCount(g_.EdgeStart(e)) == 0;
    }

    void RemoveEdge(EdgeId e) {
        VertexId start = g_.EdgeStart(e);
        VertexId end = g_.EdgeEnd(e);
        start->Lock();
        end->Lock();
        g_.RemoveEdge(e);
        start->Unlock();
        end->Unlock();
    }

public:

    ParallelTipClippingFunctor(Graph& g, size_t length_bound, HandlerF handler_f)
            : g_(g),
              length_bound_(length_bound),
              handler_f_(handler_f) {

    }

    bool operator()(VertexId v) const {
        if (LockingOutgoingCount(v) == 0)
            return false;

        vector<EdgeId> tips;
        //don't need lock here after the previous check
        for (EdgeId e : g_.IncomingEdges(v)) {
            if (IsIncomingTip(e)) {
                tips.push_back(e);
            }
        }

        //if all of edges are tips, leave the longest one
        if (tips.size() == g_.IncomingEdgeCount(v)) {
            sort(tips.begin(), tips.end(), omnigraph::LengthComparator<Graph>(g_));
            tips.pop_back();
        }

        for (EdgeId e : tips) {
            if (handler_f_) {
                handler_f_(e);
            }
            RemoveEdge(e);
        }
        v->Unlock();
        return false;
    }
};

template<class Graph>
class ParallelSimpleBRFunctor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    Graph& g_;
    size_t max_length_;
    double max_coverage_;
    double max_relative_coverage_;
    size_t max_delta_;
    double max_relative_delta_;

    bool LengthDiffCheck(size_t l1, size_t l2, size_t delta) const {
        return l1 <= l2 + delta  && l2 <= l1 + delta;
    }

    EdgeId Alternative(EdgeId e, const vector<EdgeId>& edges) const {
        size_t delta = omnigraph::CountMaxDifference(max_delta_, g_.length(e), max_relative_delta_);
        for (EdgeId candidate : boost::adaptors::reverse(edges)) {
            if (g_.EdgeEnd(candidate) == g_.EdgeEnd(e)
                    && candidate != e
                    && candidate != g_.conjugate(e)
                    && LengthDiffCheck(g_.length(candidate), g_.length(e), delta)) {
                return candidate;
            }
        }
        return EdgeId(0);
    }

    bool ProcessEdges(const vector<EdgeId>& edges) {
        for (EdgeId e : edges) {
            if (g_.length(e) <= max_length_ && math::le(g_.coverage(e), max_coverage_)) {
                EdgeId alt = Alternative(e, edges);
                if (alt != EdgeId(0)
                        && math::ge(g_.coverage(alt) * max_relative_coverage_, g_.coverage(e))) {
                    //todo is not work in multiple threads for now :)
                    //Reasons: id distribution, kmer-mapping
                    g_.Glue(e, alt);
                    return true;
                }
            }
        }
        return false;
    }

    vector<VertexId> MultiEdgeDestinations(VertexId v) const {
        vector<VertexId> answer;
        set<VertexId> destinations;
        for (EdgeId e : g_.OutgoingEdges(v)) {
            VertexId end = g_.EdgeEnd(e);
            if (destinations.count(end) > 0) {
                answer.push_back(end);
            }
            destinations.insert(end);
        }
        return answer;
    }


    VertexId SingleMultiEdgeDestination(VertexId v) const {
        vector<VertexId> dests = MultiEdgeDestinations(v);
        if (dests.size() == 1) {
            return dests.front();
        } else {
            return VertexId(0);
        }
    }

    void RemoveBulges(VertexId v) {
        bool flag = true;
        while (flag) {
            vector<EdgeId> edges(g_.out_begin(v), g_.out_end(v));
            if (edges.size() == 1)
                return;
            sort(edges.begin(), edges.end(), omnigraph::CoverageComparator<Graph>(g_));
            flag = ProcessEdges(edges);
        }
    }

    bool CheckVertex(VertexId v) const {
        bool answer;
        v->Lock();
        answer = MultiEdgeDestinations(v).size() == 1
                && MultiEdgeDestinations(g_.conjugate(v)).size() == 0;
        v->Unlock();
        return answer;
    }

    size_t MinId(VertexId v) const {
        return std::min(v.int_id(), g_.conjugate(v).int_id());
    }

    bool IsMinimal(VertexId v1, VertexId v2) const {
        return MinId(v1) < MinId(v2);
    }

public:

    ParallelSimpleBRFunctor(Graph& g,
                            size_t max_length,
                            double max_coverage,
                            double max_relative_coverage,
                            size_t max_delta,
                            double max_relative_delta,
                            boost::function<void(EdgeId)> handler_f)
            : g_(g),
              max_length_(max_length),
              max_coverage_(max_coverage),
              max_relative_coverage_(max_relative_coverage),
              max_delta_(max_delta),
              max_relative_delta_(max_relative_delta) {

    }

    bool operator()(VertexId v/*, need number of vertex for stable id distribution*/) const {
        vector<VertexId> multi_dest;
        Lock(v);
        multi_dest = MultiEdgeDestinations(v);
        Unlock(v);

        if (multi_dest.size() == 1 && IsMinimal(v, multi_dest.front())) {
            VertexId dest = multi_dest.front();
            if (CheckVertex(v) && CheckVertex(g_.conjugate(dest))) {
                Lock(v);
                Lock(dest);
                RemoveBulges(v);
                Unlock(v);
                Unlock(dest);
            }
        }
        return false;
    }

};

}

}
