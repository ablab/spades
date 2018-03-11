#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/path_extender.hpp"
namespace path_extend {

class TipSearcher {
 protected:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

 public:
    virtual boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex& vertex) const = 0;
};

class PathExtenderTipSearcher: public TipSearcher {
    using TipSearcher::ScaffoldVertex;

    const Graph& g_;
    const scaffold_graph::ScaffoldGraph& scaffold_graph_;
    shared_ptr<PathExtender> path_extender_;
    size_t edge_length_threshold_;

 public:
    PathExtenderTipSearcher(const Graph &g_,
                            const scaffold_graph::ScaffoldGraph &scaffold_graph_,
                            shared_ptr<PathExtender> path_extender_,
                            size_t edge_length_threshold_);

    boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex &vertex) const override;
};

class ScaffoldGraphGapCloser {
 protected:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
 public:
    virtual void CloseGaps(ScaffoldGraph& graph) const = 0;
};

class TipFinderGapCloser: public ScaffoldGraphGapCloser {
 protected:
    using ScaffoldGraphGapCloser::ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    shared_ptr<TipSearcher> tip_searcher_;

 public:
    TipFinderGapCloser(shared_ptr<TipSearcher> tip_searcher_);

    void CloseGaps(ScaffoldGraph &graph) const override;
};
}