//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "func/func.hpp"
#include <boost/none.hpp>
#include <atomic>
#include "assembly_graph/core/graph_iterators.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "edge_removal.hpp"
#include "func/pred.hpp"
#include "utils/logger/logger.hpp"

namespace omnigraph {

template<class Graph>
using EdgeRemovalHandlerF = std::function<void(typename Graph::EdgeId)>;

template<class Graph>
class EdgeProcessingAlgorithm {
    typedef typename Graph::EdgeId EdgeId;
    typedef func::TypedPredicate<EdgeId> ProceedConditionT;

    Graph& g_;
    bool conjugate_symmetry_;
 protected:

    Graph& g() {
        return g_;
    }

    const Graph& g() const {
        return g_;
    }

    virtual bool ProcessEdge(EdgeId e) = 0;

 public:
    EdgeProcessingAlgorithm(Graph& g,
                             bool conjugate_symmetry = false)
            : g_(g), conjugate_symmetry_(conjugate_symmetry) {

    }

    virtual ~EdgeProcessingAlgorithm() {
    }

//    bool conjugate_symmetry() const {
//        return conjugate_symmetry_;
//    }

    template<class Comparator = std::less<EdgeId>>
    bool Run(const Comparator& comp = Comparator(), ProceedConditionT proceed_condition = func::AlwaysTrue<EdgeId>()) {
        bool triggered = false;
        for (auto it = g_.SmartEdgeBegin(comp, conjugate_symmetry_); !it.IsEnd(); ++it) {
            EdgeId e = *it;
            TRACE("Current edge " << g_.str(e));
            if (!proceed_condition(e)) {
                TRACE("Stop condition was reached.");
                break;
            }

            TRACE("Processing edge " << this->g().str(e));
            triggered |= ProcessEdge(e);
        };
        return triggered;
    }

 private:
    DECL_LOGGER("EdgeProcessingAlgorithm");
};

template<class Graph>
class CountingCallback {
    typedef typename Graph::EdgeId EdgeId;
    bool report_on_destruction_;
    std::atomic<size_t> cnt_;

public:
    CountingCallback(bool report_on_destruction = false) :
            report_on_destruction_(report_on_destruction), cnt_(0) {
    }

    ~CountingCallback() {
        if (report_on_destruction_)
            Report();
    }

    void HandleDelete(EdgeId /*e*/) {
        cnt_++;
    }

    void Report() {
        TRACE(cnt_ << " edges were removed.")
        cnt_ = 0;
    }

private:
    DECL_LOGGER("CountingCallback");
};

template<class Graph>
std::function<void(typename Graph::EdgeId)> AddCountingCallback(CountingCallback<Graph>& cnt_callback, std::function<void(typename Graph::EdgeId)> handler) {
    std::function<void(typename Graph::EdgeId)> cnt_handler = std::bind(&CountingCallback<Graph>::HandleDelete, std::ref(cnt_callback), std::placeholders::_1);
    return func::CombineCallbacks<typename Graph::EdgeId>(handler, cnt_handler);
}

template<class Graph>
class EdgeRemovingAlgorithm : public EdgeProcessingAlgorithm<Graph> {
    typedef EdgeProcessingAlgorithm<Graph> base;
    typedef typename Graph::EdgeId EdgeId;

    func::TypedPredicate<EdgeId> remove_condition_;
    EdgeRemover<Graph> edge_remover_;

 protected:
    virtual bool ProcessEdge(EdgeId e) {
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
    EdgeRemovingAlgorithm(Graph& g,
                          func::TypedPredicate<EdgeId> remove_condition,
                          std::function<void (EdgeId)> removal_handler = boost::none,
                          bool conjugate_symmetry = false)
            : base(g, conjugate_symmetry),
              remove_condition_(remove_condition),
              edge_remover_(g, removal_handler) {}

 private:
    DECL_LOGGER("EdgeRemovingAlgorithm");
};

}
