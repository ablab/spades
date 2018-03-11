#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/path_extender.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {

class TipSearcher {
 protected:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

 public:
    virtual boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex& vertex) const = 0;
};

class PathExtenderTipSearcher: public TipSearcher {
    using TipSearcher::ScaffoldVertex;
    typedef unordered_map<EdgeId, ScaffoldGraphExtractor::VertexSet> edge_to_scaff_vertex_set_t;

    const Graph& g_;
    const scaffold_graph::ScaffoldGraph& scaffold_graph_;
    shared_ptr<PathExtender> path_extender_;
    edge_to_scaff_vertex_set_t edge_to_scaff_vertex_set_;
    size_t edge_length_threshold_;

 public:
    PathExtenderTipSearcher(const Graph &g_,
                            const scaffold_graph::ScaffoldGraph &scaffold_graph_,
                            shared_ptr<PathExtender> path_extender_,
                            const edge_to_scaff_vertex_set_t &edge_to_scaff_vertex_set,
                            size_t edge_length_threshold_);

    boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex &vertex) const override;

    DECL_LOGGER("PathExtenderTipSearcher");
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

    DECL_LOGGER("TipFinderGapCloser");
};
}