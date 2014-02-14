#pragma once

#include "standard_base.hpp"
#include "omni/graph_processing_algorithm.hpp"
#include "omni/basic_edge_conditions.hpp"
#include "omni/bulge_remover.hpp"
#include <boost/range/adaptor/reversed.hpp>

namespace debruijn {

namespace simplification {

//template<class VertexId>
//class LockedVertexId : boost::noncopyable {
//    typedef typename VertexId::type VertexT;
//    VertexId v_;
//
//public:
//    explicit LockedVertexId(VertexId v) : v_(v) {
//    }
//
//    VertexT *get        () const  { return v_.get();  }
//    VertexT& operator*  () const  { return *v_;     }
//    VertexT* operator-> () const  { return v_.get();}
//
//    const VertexId& vertex_id() const {
//        return v_;
//    }
//
//    bool operator==(const LockedVertexId &rhs) const {
//        return v_ == rhs.v_;
//    }
//
//    bool operator!=(const LockedVertexId &rhs) const {
//        return !operator ==(rhs);
//    }
//
//    bool operator<(const LockedVertexId &rhs) const {
//      return v_ < rhs.v_;
//    }
//
//    bool operator==(const VertexId &rhs) const {
//        return v_ == rhs;
//    }
//
//    bool operator!=(const VertexId &rhs) const {
//        return !operator ==(rhs);
//    }
//
//    bool operator<(const VertexId &rhs) const {
//      return v_ < rhs;
//    }
//
//    size_t hash() const {
//      return v_->int_id_;
//    }
//
//    size_t int_id() const {
//      return v_->int_id_;
//    }
//
//};

template<class VertexId>
class VertexLock {
    VertexId v_;
public:
    VertexLock(VertexId v) : v_(v) {
        v_->Lock();
    }

    ~VertexLock() {
        v_->Unlock();
    }
};

template<class Graph>
class ParallelTipClippingFunctor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef boost::function<void(EdgeId)> HandlerF;
    typedef VertexLock<VertexId> VertexLockT;

    Graph& g_;
    size_t length_bound_;
    HandlerF handler_f_;

    size_t LockingIncomingCount(VertexId v) const {
        VertexLockT lock(v);
        return g_.IncomingEdgeCount(v);
    }

    size_t LockingOutgoingCount(VertexId v) const {
        VertexLockT lock(v);
        return g_.OutgoingEdgeCount(v);
    }

    bool IsIncomingTip(EdgeId e) const {
        return g_.length(e) <= length_bound_ 
                && LockingIncomingCount(g_.EdgeStart(e)) + LockingOutgoingCount(g_.EdgeStart(e)) == 1;
    }

    void RemoveEdge(EdgeId e) {
        VertexLockT lock1(g_.EdgeStart(e));
        VertexLockT lock2(g_.EdgeEnd(e));
        g_.DeleteEdge(e);
    }

public:

    ParallelTipClippingFunctor(Graph& g, size_t length_bound, HandlerF handler_f = 0)
            : g_(g),
              length_bound_(length_bound),
              handler_f_(handler_f) {

    }

    bool operator()(VertexId v) {
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
        if (!tips.empty() && tips.size() == g_.IncomingEdgeCount(v)) {
            sort(tips.begin(), tips.end(), omnigraph::LengthComparator<Graph>(g_));
            tips.pop_back();
        }

        for (EdgeId e : tips) {
            if (handler_f_) {
                handler_f_(e);
            }
            RemoveEdge(e);
        }
        return false;
    }
};

template<class Graph>
class ParallelSimpleBRFunctor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef VertexLock<VertexId> VertexLockT;

    Graph& g_;
    size_t max_length_;
    double max_coverage_;
    double max_relative_coverage_;
    size_t max_delta_;
    double max_relative_delta_;
    boost::function<void(EdgeId)> handler_f_;

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
                    handler_f_(e);
                    g_.GlueEdges(e, alt);
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
        VertexLockT lock(v);
        return MultiEdgeDestinations(v).size() == 1
                && MultiEdgeDestinations(g_.conjugate(v)).size() == 0;
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
                            boost::function<void(EdgeId)> handler_f = 0)
            : g_(g),
              max_length_(max_length),
              max_coverage_(max_coverage),
              max_relative_coverage_(max_relative_coverage),
              max_delta_(max_delta),
              max_relative_delta_(max_relative_delta),
              handler_f_(handler_f) {

    }

    bool operator()(VertexId v/*, need number of vertex for stable id distribution*/) {
        vector<VertexId> multi_dest;

        {
            VertexLockT lock(v);
            multi_dest = MultiEdgeDestinations(v);
        }

        if (multi_dest.size() == 1 && IsMinimal(v, multi_dest.front())) {
            VertexId dest = multi_dest.front();
            if (CheckVertex(v) && CheckVertex(g_.conjugate(dest))) {
                VertexLockT lock1(v);
                VertexLockT lock2(dest);
                RemoveBulges(v);
            }
        }
        return false;
    }

};

//currently just a stab and not parallel at all,
//no way to write it parallel before deletion of edge requires two locks
template<class Graph>
class ParallelLowCoverageFunctor {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef VertexLock<VertexId> VertexLockT;
    typedef boost::function<void(EdgeId)> HandlerF;

    Graph& g_;
    shared_ptr<func::Predicate<EdgeId>> ec_condition_;
    HandlerF handler_f_;

    vector<EdgeId> edges_to_remove_;

public:

    ParallelLowCoverageFunctor(Graph& g,
                               size_t max_length,
                               double max_coverage,
                               HandlerF handler_f = 0)
            : g_(g),
              ec_condition_(
                      func::And<EdgeId>(
                      func::And<EdgeId>(make_shared<omnigraph::LengthUpperBound<Graph>>(g, max_length),
                                        make_shared<omnigraph::CoverageUpperBound<Graph>>(g, max_coverage)),
                      make_shared<omnigraph::AlternativesPresenceCondition<Graph>>(g))),
              handler_f_(handler_f)
    {

    }

    bool operator()(EdgeId e) {
        if (ec_condition_->Check(e)) {
            edges_to_remove_.push_back(e);
        }
        return false;
    }

    void RemoveCollectedEdges() {
        omnigraph::SmartSetIterator<Graph, EdgeId> to_delete(g_, edges_to_remove_.begin(), edges_to_remove_.end());
        while (!to_delete.IsEnd()) {
            EdgeId e = *to_delete;
            handler_f_(e);
            g_.DeleteEdge(e);
            ++to_delete;
        }
    }

};

}

}
