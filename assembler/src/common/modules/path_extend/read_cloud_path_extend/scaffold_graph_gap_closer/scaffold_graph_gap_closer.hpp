#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/path_extender.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace path_extend {

class TipSearcher {
 protected:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

 public:
    virtual boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex &vertex) const = 0;
};

class TipExtender {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    const Graph &g_;
    shared_ptr<PathExtender> path_extender_;
    size_t edge_length_threshold_;

 public:
    TipExtender(const Graph &g, shared_ptr<PathExtender> path_extender, size_t edge_length_threshold);

    BidirectionalPath ExtendTip(const ScaffoldVertex &vertex) const;
};

class PathExtenderTipSearcher : public TipSearcher {
    using TipSearcher::ScaffoldVertex;
    typedef unordered_map<EdgeId, ScaffoldGraphExtractor::VertexSet> edge_to_scaff_vertex_set_t;

    const Graph &g_;
    shared_ptr<TipExtender> tip_extender_;
    edge_to_scaff_vertex_set_t edge_to_scaff_vertex_set_;
    size_t edge_length_threshold_;

 public:
    PathExtenderTipSearcher(const Graph &g_,
                            shared_ptr<TipExtender> path_extender_,
                            const edge_to_scaff_vertex_set_t &edge_to_scaff_vertex_set,
                            size_t edge_length_threshold_);

    boost::optional<ScaffoldVertex> FindTip(const ScaffoldVertex &vertex) const override;

    DECL_LOGGER("PathExtenderTipSearcher");
};

class ScaffoldGraphGapCloser {
 protected:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
 public:
    virtual void CloseGaps(ScaffoldGraph &graph) const = 0;
};

class TipFinderGapCloser : public ScaffoldGraphGapCloser {
 protected:
    using ScaffoldGraphGapCloser::ScaffoldGraph;
    using ScaffoldGraphGapCloser::ScaffoldVertex;
    shared_ptr<TipSearcher> tip_searcher_;

 public:
    TipFinderGapCloser(shared_ptr<TipSearcher> tip_searcher_);

    void CloseGaps(ScaffoldGraph &graph) const override;

    DECL_LOGGER("TipFinderGapCloser");
};

struct TipPair {
  typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

  ScaffoldVertex first_;
  ScaffoldVertex second_;
  const double score_;

  TipPair(const ScaffoldVertex &first_, const ScaffoldVertex &second_, const double score_);
};

class ScoreFunctionGapCloser : public ScaffoldGraphGapCloser {
 protected:
    using ScaffoldGraphGapCloser::ScaffoldGraph;
    using ScaffoldGraphGapCloser::ScaffoldVertex;
    typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    const Graph &g_;
    shared_ptr<TipExtender> tip_extender_;
    shared_ptr<BarcodeEntryCollector> entry_collector_;
    const double score_threshold_;

 public:
    ScoreFunctionGapCloser(const Graph &g_,
                           shared_ptr<TipExtender> tip_extender_,
                           shared_ptr<BarcodeEntryCollector> entry_collector_,
                           double score_threshold_);

    void CloseGaps(ScaffoldGraph &graph) const override;

    DECL_LOGGER("ScoreFunctionGapCloser");
};
}