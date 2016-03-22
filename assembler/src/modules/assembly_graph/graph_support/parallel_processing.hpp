//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "dev_support/logger/logger.hpp"
#include "assembly_graph/graph_core/graph_iterators.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "dev_support/openmp_wrapper.h"

namespace omnigraph {

template<class ItVec, class SmartIt, class Predicate>
void FillInterestingFromChunkIterators(const ItVec& chunk_iterators,
                                       SmartIt& smart_it,
                                       const Predicate& predicate) {
    VERIFY(chunk_iterators.size() > 1);
    typedef typename Predicate::checked_type ElementType;
    std::vector<std::vector<ElementType>> of_interest(omp_get_max_threads());

    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
        for (auto it = chunk_iterators[i], end = chunk_iterators[i + 1]; it != end; ++it) {
            ElementType t = *it;
            if (predicate(t)) {
                of_interest[omp_get_thread_num()].push_back(t);
            }
        }
    }

    for (auto& chunk : of_interest) {
        smart_it.insert(chunk.begin(), chunk.end());
        chunk.clear();
    }
}

template<class Graph, class ElementId = typename Graph::EdgeId>
class TrivialInterestingElementFinder {
public:

    TrivialInterestingElementFinder() {
    }

    template<class SmartIt>
    bool Run(SmartIt& /*it*/) const {
        return false;
    }
};

template<class Graph, class ElementId = typename Graph::EdgeId>
class SimpleInterestingElementFinder {
    typedef GraphEdgeIterator<Graph> EdgeIt;

    const Graph& g_;
    pred::TypedPredicate<ElementId> condition_;
public:

    SimpleInterestingElementFinder(const Graph& g,
                                   pred::TypedPredicate<ElementId> condition = pred::AlwaysTrue<ElementId>())
            :  g_(g), condition_(condition) {}

    template<class SmartIt>
    bool Run(SmartIt& interest) const {
        for (EdgeIt it = EdgeIt(g_, g_.begin()), end = EdgeIt(g_, g_.end()); it != end; ++it) {
            if (condition_(*it)) {
                interest.push(*it);
            }
        }
        return false;
    }
};

template<class Graph, class ElementId = typename Graph::EdgeId>
class ParallelInterestingElementFinder {
    typedef GraphEdgeIterator<Graph> EdgeIt;

    const Graph& g_;
    pred::TypedPredicate<ElementId> condition_;
    const size_t chunk_cnt_;
public:

    ParallelInterestingElementFinder(const Graph& g,
                                     pred::TypedPredicate<ElementId> condition,
                                     size_t chunk_cnt)
            : g_(g), condition_(condition), chunk_cnt_(chunk_cnt) {}

    template<class SmartIt>
    bool Run(SmartIt& it) const {
        TRACE("Looking for interesting elements");
        TRACE("Splitting graph into " << chunk_cnt_ << " chunks");
        FillInterestingFromChunkIterators(IterationHelper<Graph, ElementId>(g_).Chunks(chunk_cnt_), it, condition_);
        TRACE("Found " << it.size() << " interesting elements");
        return false;
    }
private:
    DECL_LOGGER("ParallelInterestingElementFinder");
};

template<class Graph>
class PersistentAlgorithmBase {
    Graph& g_;
protected:

    PersistentAlgorithmBase(Graph& g) : g_(g) {}

    Graph& g() { return g_; }
    const Graph& g() const { return g_; }
public:
    virtual ~PersistentAlgorithmBase() {}
    virtual bool Run(bool force_primary_launch = false) = 0;
};

//todo use add_condition in it_
template<class Graph, class ElementId, class InterestingElementFinder,
         class Comparator = std::less<ElementId>>
class PersistentProcessingAlgorithm : public PersistentAlgorithmBase<Graph> {
    InterestingElementFinder interest_el_finder_;

    SmartSetIterator<Graph, ElementId, Comparator> it_;
    //todo remove
    bool tracking_;
    size_t total_iteration_estimate_;

    size_t curr_iteration_;

protected:

    virtual bool Process(ElementId el) = 0;
    virtual bool Proceed(ElementId /*el*/) const { return true; }

    virtual void PrepareIteration(size_t /*it_cnt*/, size_t /*total_it_estimate*/) {}

public:

    PersistentProcessingAlgorithm(Graph& g,
                                      const InterestingElementFinder& interest_el_finder,
                                      bool canonical_only = false,
                                      const Comparator& comp = Comparator(),
                                      bool track_changes = true,
                                      size_t total_iteration_estimate = -1ul) :
                                      PersistentAlgorithmBase<Graph>(g),
                                      interest_el_finder_(interest_el_finder),
                                      it_(g, true, comp, canonical_only),
                                      tracking_(track_changes),
                                      total_iteration_estimate_(total_iteration_estimate),
                                      curr_iteration_(0) {
        it_.Detach();
    }

    bool Run(bool force_primary_launch = false) {
        bool primary_launch = !tracking_ || (curr_iteration_ == 0) || force_primary_launch ;
        if (!it_.IsAttached()) {
            it_.Attach();
        }
        if (primary_launch) {
            it_.clear();
            TRACE("Primary launch.");
            TRACE("Start preprocessing");
            interest_el_finder_.Run(it_);
            TRACE(it_.size() << " edges to process after preprocessing");
        } else {
            TRACE(it_.size() << " edges to process");
            VERIFY(tracking_);
        }

        if (curr_iteration_ >= total_iteration_estimate_) {
            PrepareIteration(total_iteration_estimate_ - 1, total_iteration_estimate_);
        } else {
            PrepareIteration(curr_iteration_, total_iteration_estimate_);
        }

        bool triggered = false;
        TRACE("Start processing");
        for (; !it_.IsEnd(); ++it_) {
            ElementId el = *it_;
            if (!Proceed(el)) {
                TRACE("Proceed condition turned false on element " << this->g().str(el));
                it_.ReleaseCurrent();
                break;
            }
            TRACE("Processing edge " << this->g().str(el));
            triggered |= Process(el);
        }
        TRACE("Finished processing. Triggered = " << triggered);
        if (!tracking_)
            it_.Detach();

        curr_iteration_++;
        return triggered;
    }

};

template<class Graph, class InterestingEdgeFinder,
         class Comparator = std::less<typename Graph::EdgeId>>
class PersistentEdgeRemovingAlgorithm : public PersistentProcessingAlgorithm<Graph,
                                                                            typename Graph::EdgeId,
                                                                            InterestingEdgeFinder, Comparator> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentProcessingAlgorithm<Graph, EdgeId, InterestingEdgeFinder, Comparator> base;
    EdgeRemover<Graph> edge_remover_;
public:
    PersistentEdgeRemovingAlgorithm(Graph& g,
                                    const InterestingEdgeFinder& interest_edge_finder,
                                    std::function<void(EdgeId)> removal_handler = boost::none,
                                    bool canonical_only = false,
                                    const Comparator& comp = Comparator(),
                                    bool track_changes = true,
                                    size_t total_iteration_estimate = -1ul)
            : base(g, interest_edge_finder,
                   canonical_only, comp, track_changes,
                   total_iteration_estimate),
                   edge_remover_(g, removal_handler) {

    }

protected:

    virtual bool ShouldRemove(EdgeId e) const = 0;

    bool Process(EdgeId e) override {
        TRACE("Checking edge " << this->g().str(e) << " for the removal condition");
        if (ShouldRemove(e)) {
            TRACE("Check passed, removing");
            edge_remover_.DeleteEdge(e);
            return true;
        }
        TRACE("Check not passed");
        return false;
    }

};

template<class Graph, class InterestingEdgeFinder,
         class Comparator = std::less<typename Graph::EdgeId>>
class ConditionEdgeRemovingAlgorithm : public PersistentEdgeRemovingAlgorithm<Graph,
                                                                              InterestingEdgeFinder, Comparator> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentEdgeRemovingAlgorithm<Graph, InterestingEdgeFinder, Comparator> base;
    pred::TypedPredicate<EdgeId> remove_condition_;
protected:

    bool ShouldRemove(EdgeId e) const override {
        return remove_condition_(e);
    }

public:
    ConditionEdgeRemovingAlgorithm(Graph& g,
                                   const InterestingEdgeFinder& interest_edge_finder,
                                   pred::TypedPredicate<EdgeId> remove_condition,
                                   std::function<void(EdgeId)> removal_handler = boost::none,
                                   bool canonical_only = false,
                                   const Comparator& comp = Comparator(),
                                   bool track_changes = true)
            : base(g, interest_edge_finder,
                   removal_handler,
                   canonical_only, comp, track_changes),
                   remove_condition_(remove_condition) {

    }
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class ParallelEdgeRemovingAlgorithm : public ConditionEdgeRemovingAlgorithm<Graph,
                                                ParallelInterestingElementFinder<Graph>, Comparator> {
    typedef ConditionEdgeRemovingAlgorithm<Graph,
            ParallelInterestingElementFinder<Graph>, Comparator> base;
    typedef typename Graph::EdgeId EdgeId;

public:
    ParallelEdgeRemovingAlgorithm(Graph& g,
                                  pred::TypedPredicate<EdgeId> remove_condition,
                                  size_t chunk_cnt,
                                  std::function<void(EdgeId)> removal_handler = boost::none,
                                  bool canonical_only = false,
                                  const Comparator& comp = Comparator(),
                                  bool track_changes = true)
            : base(g,
                   ParallelInterestingElementFinder<Graph>(g, remove_condition, chunk_cnt),
                   remove_condition, removal_handler,
                   canonical_only, comp, track_changes) {
    }

};

}
