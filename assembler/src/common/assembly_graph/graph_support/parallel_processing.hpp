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
#include "utils/parallel/openmp_wrapper.h"

namespace omnigraph {

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

template<class Graph, class ElementId>
using InterestingFinderPtr = std::shared_ptr<InterestingElementFinder<Graph, ElementId>>;

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

    template<class ItVec, class Condition, class Handler>
    static void FindInterestingFromChunkIterators(const ItVec& chunk_iterators,
                                           const Condition& predicate,
                                           const Handler& handler) {
        VERIFY(chunk_iterators.size() > 1);
        DEBUG("Parallel search for elements of interest");
        typedef typename Condition::checked_type ElementType;
        std::vector<std::vector<ElementType>> of_interest(chunk_iterators.size() - 1);

        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < chunk_iterators.size() - 1; ++i) {
            DEBUG("Processing chunk " << i << " by thread " << omp_get_thread_num());
            size_t cnt = 0;
            for (auto it = chunk_iterators[i], end = chunk_iterators[i + 1]; it != end; ++it) {
                ElementType t = *it;
                if (predicate(t)) {
                    of_interest[i].push_back(t);
                }
                cnt++;
            }
            DEBUG("Processed chunk " << i << ". " << cnt << " elements identified as potential candidates");
        }

        DEBUG("Merging chunks");
        for (auto& chunk : of_interest) {
            for (const auto& el : chunk) {
                handler(el);
            }
            chunk.clear();
        }
        DEBUG("Chunks merged");
    }

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

    /**
     * Launches graph processing
     * @param force_primary_launch flag forcing the refilling of the set of elements of interest
     * @param iter_run_progress progress coefficient for iterative algorithms passed here
     * @return number of trigger events
     */
    virtual size_t Run(bool force_primary_launch = false,
                       double iter_run_progress = 1.) = 0;

protected:
    DECL_LOGGER("Simplification");
};

template<class Graph>
using AlgoPtr = std::shared_ptr<omnigraph::PersistentAlgorithmBase<Graph>>;

template<class Graph>
using EdgeConditionT = func::TypedPredicate<typename Graph::EdgeId>;

template<class Graph>
class AlgorithmRunningHelper {
    typedef PersistentAlgorithmBase<Graph> Algo;
public:
    static size_t RunAlgo(Algo &algo, const string &comment = "",
                 bool force_primary_launch = false,
                 double iter_run_progress = 1.) {
        if (!comment.empty()) {INFO("Running " << comment);}
        size_t triggered = algo.Run(force_primary_launch, iter_run_progress);
        if (!comment.empty()) {INFO(comment << " triggered " << triggered << " times");}
        return triggered;
    }

    static size_t RunAlgo(AlgoPtr<Graph> algo_ptr, const string &comment = "",
                          bool force_primary_launch = false,
                          double iter_run_progress = 1.) {
        if (algo_ptr)
            return RunAlgo(*algo_ptr, comment, force_primary_launch, iter_run_progress);
        return 0;
    }

    static size_t IterativeThresholdsRun(Algo &algo,
                                         const size_t iteration_cnt = 1,
                                         bool all_primary = false,
                                         bool first_primary = true) {
        size_t total_triggered = 0;
        for (size_t i = 0; i < iteration_cnt; ++i) {
            DEBUG("Iteration " << i);
            size_t algo_triggered = algo.Run(all_primary || (i == 0 && first_primary),
                                double(i + 1) / double(iteration_cnt));
            DEBUG("Triggered " << algo_triggered << " times on iteration " << (i + 1));
            total_triggered += algo_triggered;
        }
        return total_triggered;
    }

    //TODO use enum instead all_primary/first_primary flags
    static size_t LoopedRun(Algo &algo, size_t min_it_cnt = 1,
                            size_t max_it_cnt = size_t(-1),
                            bool all_primary = false,
                            bool first_primary = true,
                            double iter_run_progress = 1.) {
        size_t triggered = 0;

        bool changed = true;
        for (size_t i = 0; i < max_it_cnt && (changed || i < min_it_cnt); ++i) {
            bool primary = all_primary || (first_primary && i == 0);
            DEBUG("Iteration " << (i + 1));
            size_t algo_triggered = algo.Run(primary, iter_run_progress);
            DEBUG("Triggered " << algo_triggered << " times on iteration " << (i + 1));
            changed = (algo_triggered > 0);
            triggered += algo_triggered;
        }
        return triggered;
    }

    static size_t LoopedRunPrimaryOpening(Algo &algo, size_t primary_cnt,
                                          size_t max_it_cnt = size_t(-1),
                                          double iter_run_progress = 1.) {
        VERIFY(primary_cnt > 0 && max_it_cnt >= primary_cnt);
        size_t triggered = 0;
        if (primary_cnt > 1)
            triggered += LoopedRun(algo, primary_cnt - 1, primary_cnt - 1, true, true, iter_run_progress);
        //parameter defaults used to get to the last one
        triggered += LoopedRun(algo, 1, max_it_cnt - primary_cnt + 1, false, true, iter_run_progress);
        return triggered;
    }

private:
    DECL_LOGGER("Simplification");
};

template<class Graph>
class AdapterAlgorithm : public PersistentAlgorithmBase<Graph> {
    std::function<size_t ()> func_;

public:
    template <class F>
    AdapterAlgorithm(Graph& g, F f) :
            PersistentAlgorithmBase<Graph>(g) {
        func_ = [=] () {return (size_t) f();};
    }

    size_t Run(bool, double) override {
        return func_();
    }
};

template<class Graph>
class CompositeAlgorithm : public PersistentAlgorithmBase<Graph> {
    std::vector<std::pair<AlgoPtr<Graph>, std::string>> algos_;
    //TODO consider moving up hierarchy or some other solution
    std::function<void ()> launch_callback_;

public:
    CompositeAlgorithm(Graph& g,
                       std::function<void ()> launch_callback = nullptr) :
            PersistentAlgorithmBase<Graph>(g),
            launch_callback_(launch_callback) {
    }

    template<typename Algo, typename... Args>
    void AddAlgo(const std::string &desc, Args&&... args) {
        AddAlgo(std::make_shared<Algo>(std::forward<Args>(args)...), desc);
    };

    void AddAlgo(AlgoPtr<Graph> algo, const std::string &desc = "") {
        if (algo)
            algos_.push_back(std::make_pair(algo, desc));
    }

    size_t Run(bool force_primary_launch = false,
               double iter_run_progress = 1.) override {
        if (launch_callback_)
            launch_callback_();
        size_t triggered = 0;

        for (const auto &algo_info : algos_) {
            triggered += AlgorithmRunningHelper<Graph>::RunAlgo(*algo_info.first,
                                                                algo_info.second,
                                                                force_primary_launch,
                                                                iter_run_progress);
        }
        return triggered;
    }

};

template<class Graph>
class LoopedAlgorithm : public PersistentAlgorithmBase<Graph> {
    AlgoPtr<Graph> algo_;
    size_t min_iter_cnt_;
    size_t max_iter_cnt_;
    bool force_primary_for_all_;

public:
    LoopedAlgorithm(Graph& g, AlgoPtr<Graph> algo,
                    size_t min_iter_cnt = 1,
                    size_t max_iter_cnt = size_t(-1),
                    bool force_primary_for_all = false) :
            PersistentAlgorithmBase<Graph>(g),
            algo_(algo),
            min_iter_cnt_(min_iter_cnt),
            max_iter_cnt_(max_iter_cnt),
            force_primary_for_all_(force_primary_for_all) {
        VERIFY(min_iter_cnt > 0);
    }

    size_t Run(bool force_primary_launch = false,
               double iter_run_progress = 1.) override {
        return AlgorithmRunningHelper<Graph>::LoopedRun(*algo_, min_iter_cnt_, max_iter_cnt_,
                                                        force_primary_launch && force_primary_for_all_,
                                                        force_primary_launch, iter_run_progress);
    }

};

//FIXME only potentially relevant edges should be stored at any point
template<class Graph, class ElementId,
         class Comparator = std::less<ElementId>>
class PersistentProcessingAlgorithm : public PersistentAlgorithmBase<Graph> {
protected:
    typedef std::shared_ptr<InterestingElementFinder<Graph, ElementId>> CandidateFinderPtr;
    CandidateFinderPtr interest_el_finder_;

private:
    SmartSetIterator<Graph, ElementId, Comparator> it_;
    const bool tracking_;

protected:
    void ReturnForConsideration(ElementId el) {
        it_.push(el);
    }

    virtual bool Process(ElementId el) = 0;
    virtual bool Proceed(ElementId /*el*/) const { return true; }
    virtual void PrepareIteration(double /*iter_run_progress*/ = 1.) {}

public:

    PersistentProcessingAlgorithm(Graph& g,
                                  CandidateFinderPtr interest_el_finder,
                                  bool canonical_only = false,
                                  const Comparator& comp = Comparator(),
                                  bool track_changes = true) :
            PersistentAlgorithmBase<Graph>(g),
            interest_el_finder_(interest_el_finder),
            it_(g, true, comp, canonical_only),
            tracking_(track_changes) {
        it_.Detach();
    }

    size_t Run(bool force_primary_launch = false,
               double iter_run_progress = 1.) override {
        bool primary_launch = force_primary_launch ;
        if (!it_.IsAttached()) {
            it_.Attach();
            primary_launch = true;
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

        //PrepareIteration(std::min(curr_iteration_, total_iteration_estimate_ - 1), total_iteration_estimate_);
        PrepareIteration(iter_run_progress);

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

//TODO use coverage order?
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
              //condition_(second_check ? condition : func::AlwaysTrue<EdgeId>()),
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
