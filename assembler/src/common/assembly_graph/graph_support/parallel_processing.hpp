//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/logger/logger.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include "assembly_graph/graph_support/graph_processing_algorithm.hpp"
#include "utils/openmp_wrapper.h"

namespace omnigraph {

template<class ItVec, class Condition, class Handler>
void FindInterestingFromChunkIterators(const ItVec& chunk_iterators,
                                       const Condition& predicate,
                                       const Handler& handler) {
    VERIFY(chunk_iterators.size() > 1);
    typedef typename Condition::checked_type ElementType;
    std::vector<std::vector<ElementType>> of_interest(omp_get_max_threads());

    #pragma omp parallel for schedule(guided)
    for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
        size_t cnt = 0;
        for (auto it = chunk_iterators[i], end = chunk_iterators[i + 1]; it != end; ++it) {
             ElementType t = *it;
             if (predicate(t)) {
                 of_interest[omp_get_thread_num()].push_back(t);
             }
             cnt++;
         }
         DEBUG("Processed " << cnt << " elements as potential candidates by thread " << omp_get_thread_num());
    }

    for (auto& chunk : of_interest) {
        for (const auto& el : chunk) {
            handler(el);
        }
        chunk.clear();
    }
}

template<class Graph, class ElementId>
class InterestingElementFinder {
protected:
    typedef std::function<void (ElementId)> HandlerF;
    const func::TypedPredicate<ElementId> condition_;
public:

    InterestingElementFinder(func::TypedPredicate<ElementId> condition):
            condition_(condition) {
    }

    virtual ~InterestingElementFinder() {}

    virtual bool Run(const Graph& /*g*/, HandlerF /*handler*/) const = 0;
};

template<class Graph, class ElementId = typename Graph::EdgeId>
class TrivialInterestingElementFinder :
        public InterestingElementFinder<Graph, ElementId> {
public:

    TrivialInterestingElementFinder() :
            InterestingElementFinder<Graph, ElementId>(func::AlwaysTrue<ElementId>()) {
    }

    bool Run(const Graph& /*g*/, std::function<void (ElementId)> /*handler*/) const override {
        return false;
    }
};

template<class Graph, class ElementId = typename Graph::EdgeId>
class SimpleInterestingElementFinder : public InterestingElementFinder<Graph, ElementId> {
    typedef InterestingElementFinder<Graph, ElementId> base;
    typedef typename base::HandlerF HandlerF;
public:

    SimpleInterestingElementFinder(func::TypedPredicate<ElementId> condition = func::AlwaysTrue<ElementId>())
            :  base(condition) {}

    bool Run(const Graph& g, HandlerF handler) const override {
        const IterationHelper<Graph, ElementId> it_helper(g);
        for (auto it = it_helper.begin(), end = it_helper.end(); it != end; ++it) {
            if (this->condition_(*it)) {
                handler(*it);
            }
        }
        return false;
    }
};

template<class Graph, class ElementId = typename Graph::EdgeId>
class ParallelInterestingElementFinder : public InterestingElementFinder<Graph, ElementId> {
    typedef InterestingElementFinder<Graph, ElementId> base;
    typedef typename base::HandlerF HandlerF;

    const size_t chunk_cnt_;
public:

    ParallelInterestingElementFinder(func::TypedPredicate<ElementId> condition,
                                     size_t chunk_cnt)
            : base(condition), chunk_cnt_(chunk_cnt) {}

    bool Run(const Graph& g, HandlerF handler) const override {
        TRACE("Looking for interesting elements");
        TRACE("Splitting graph into " << chunk_cnt_ << " chunks");
        FindInterestingFromChunkIterators(IterationHelper<Graph, ElementId>(g).Chunks(chunk_cnt_),
                                          this->condition_, handler);
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
    virtual size_t Run(bool force_primary_launch = false) = 0;
};

template<class Algo>
inline size_t LoopedRun(Algo& algo) {
    size_t total_triggered = 0;
    bool run = true;
    while (run) {
        size_t triggered = algo.Run();
        total_triggered += triggered;
        run = (triggered > 0);
    }
    return total_triggered;
}

//todo only potentially relevant edges should be stored at any point
template<class Graph, class ElementId,
         class Comparator = std::less<ElementId>>
class PersistentProcessingAlgorithm : public PersistentAlgorithmBase<Graph> {
protected:
    typedef std::shared_ptr<InterestingElementFinder<Graph, ElementId>> CandidateFinderPtr;
    CandidateFinderPtr interest_el_finder_;

private:
    SmartSetIterator<Graph, ElementId, Comparator> it_;
    bool tracking_;
    size_t total_iteration_estimate_;
    size_t curr_iteration_;

protected:
    void ReturnForConsideration(ElementId el) {
        it_.push(el);
    }

    virtual bool Process(ElementId el) = 0;
    virtual bool Proceed(ElementId /*el*/) const { return true; }

    virtual void PrepareIteration(size_t /*it_cnt*/, size_t /*total_it_estimate*/) {}

public:

    PersistentProcessingAlgorithm(Graph& g,
                                  const CandidateFinderPtr& interest_el_finder,
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

    size_t Run(bool force_primary_launch = false) override {
        bool primary_launch = !tracking_ || (curr_iteration_ == 0) || force_primary_launch ;
        if (!it_.IsAttached()) {
            it_.Attach();
        }
        if (primary_launch) {
            it_.clear();
            TRACE("Primary launch.");
            TRACE("Start searching for relevant elements");
            interest_el_finder_->Run(this->g(), [&](ElementId el) {it_.push(el);});
            TRACE(it_.size() << " elements to consider");
        } else {
            TRACE(it_.size() << " elements to consider");
            VERIFY(tracking_);
        }

        PrepareIteration(std::min(curr_iteration_, total_iteration_estimate_ - 1), total_iteration_estimate_);

        size_t triggered = 0;
        TRACE("Start processing");
        for (; !it_.IsEnd(); ++it_) {
            ElementId el = *it_;
            if (!Proceed(el)) {
                TRACE("Proceed condition turned false on element " << this->g().str(el));
                it_.ReleaseCurrent();
                break;
            }
            TRACE("Processing edge " << this->g().str(el));
            if (Process(el))
                triggered++;
        }
        TRACE("Finished processing. Triggered = " << triggered);
        if (!tracking_)
            it_.Detach();

        curr_iteration_++;
        return triggered;
    }
private:
    DECL_LOGGER("PersistentProcessingAlgorithm"); 
};

template<class Graph,
        class Comparator = std::less<typename Graph::EdgeId>>
class ParallelEdgeRemovingAlgorithm : public PersistentProcessingAlgorithm<Graph,
        typename Graph::EdgeId,
        Comparator> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentProcessingAlgorithm<Graph, EdgeId, Comparator> base;

    const func::TypedPredicate<EdgeId> remove_condition_;
    EdgeRemover<Graph> edge_remover_;

protected:

    bool Process(EdgeId e) override {
        TRACE("Checking edge " << this->g().str(e) << " for the removal condition");
        if (remove_condition_(e)) {
            TRACE("Check passed, removing");
            edge_remover_.DeleteEdge(e);
            return true;
        }
        TRACE("Check not passed");
        return false;
    }

public:
    ParallelEdgeRemovingAlgorithm(Graph& g,
                                  func::TypedPredicate<EdgeId> remove_condition,
                                  size_t chunk_cnt,
                                  std::function<void(EdgeId)> removal_handler = boost::none,
                                  bool canonical_only = false,
                                  const Comparator& comp = Comparator(),
                                  bool track_changes = true)
            : base(g,
                   std::make_shared<ParallelInterestingElementFinder<Graph>>(remove_condition, chunk_cnt),
                   canonical_only, comp, track_changes),
                   remove_condition_(remove_condition),
                   edge_remover_(g, removal_handler) {
    }

private:
    DECL_LOGGER("ParallelEdgeRemovingAlgorithm");
};

template<class Graph, class Comparator = std::less<typename Graph::EdgeId>>
class DisconnectionAlgorithm : public PersistentProcessingAlgorithm<Graph,
        typename Graph::EdgeId,
        Comparator> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentProcessingAlgorithm<Graph, EdgeId, Comparator> base;
    func::TypedPredicate<EdgeId> condition_;
    EdgeDisconnector<Graph> disconnector_;

public:
    DisconnectionAlgorithm(Graph& g,
                           func::TypedPredicate<EdgeId> condition,
                           size_t chunk_cnt,
                           EdgeRemovalHandlerF<Graph> removal_handler,
                           const Comparator& comp = Comparator(),
                           bool track_changes = true)
            : base(g,
                   std::make_shared<omnigraph::ParallelInterestingElementFinder<Graph>>(condition, chunk_cnt),
            /*canonical_only*/false, comp, track_changes),
              condition_(condition),
              disconnector_(g, removal_handler) {
    }

    bool Process(EdgeId e) override {
        if (condition_(e)) {
            disconnector_(e);
            return true;
        }
        return false;
    }

};


}
