//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "barcode_index/scaffold_vertex_index.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/pipeline/launch_support.hpp"
#include "modules/path_extend/extension_chooser.hpp"
#include "modules/path_extend/scaffolder2015/connection_condition2015.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldEdgePredicate : public func::AbstractPredicate<const scaffold_graph::ScaffoldGraph::ScaffoldEdge &> {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

    virtual ~ScaffoldEdgePredicate() = default;
};

class SimpleSearcher {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

    struct VertexWithDistance {
      ScaffoldVertex vertex;
      size_t distance;
      VertexWithDistance(const ScaffoldVertex &vertex, size_t distance);
    };

    SimpleSearcher(const scaffold_graph::ScaffoldGraph &graph_, const Graph &g, size_t distance_);

    std::vector<ScaffoldVertex> GetReachableVertices(const ScaffoldVertex &vertex,
                                                     const ScaffoldGraph::ScaffoldEdge &restricted_edge);
    void ProcessVertex(std::queue<VertexWithDistance> &vertex_queue, const VertexWithDistance &vertex,
                       std::unordered_set<ScaffoldVertex> &visited, const ScaffoldGraph::ScaffoldEdge &restricted_edge);
    bool AreEqual(const ScaffoldGraph::ScaffoldEdge &first, const ScaffoldGraph::ScaffoldEdge &second);

  private:
    const ScaffoldGraph &scaff_graph_;
    const Graph &g_;
    size_t distance_threshold_;

    DECL_LOGGER("SimpleSearcher");
};

class TransitiveEdgesPredicate : public ScaffoldEdgePredicate {
  public:
    using ScaffoldEdgePredicate::ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    TransitiveEdgesPredicate(const scaffold_graph::ScaffoldGraph &graph, const Graph &g, size_t distance_threshold);

    bool Check(const ScaffoldEdge &scaffold_edge) const override;

  private:
    const scaffold_graph::ScaffoldGraph scaffold_graph_;
    const Graph &g_;
    size_t distance_threshold_;

    DECL_LOGGER("TransitiveEdgesPredicate");
};

class ScaffoldEdgeScoreFunction {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    virtual double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const = 0;
    virtual ~ScaffoldEdgeScoreFunction() = default;
};

class AbstractBarcodeScoreFunction : public ScaffoldEdgeScoreFunction {
  public:
    AbstractBarcodeScoreFunction(
        const Graph &graph_,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_);

  protected:
    const Graph &graph_;
    std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;
};

class NormalizedBarcodeScoreFunction : public AbstractBarcodeScoreFunction {
  public:
    NormalizedBarcodeScoreFunction(const Graph &graph_,
                                   std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_);

    double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const override;

  protected:
    using AbstractBarcodeScoreFunction::barcode_extractor_;
    using AbstractBarcodeScoreFunction::graph_;

    DECL_LOGGER("NormalizedBarcodeScoreFunction");
};

class TrivialBarcodeScoreFunction : public AbstractBarcodeScoreFunction {
  public:
    TrivialBarcodeScoreFunction(
        const Graph &graph_,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
        size_t read_count_threshold_,
        size_t tail_threshold_);

    double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const override;

  protected:
    using AbstractBarcodeScoreFunction::barcode_extractor_;
    using AbstractBarcodeScoreFunction::graph_;
    const size_t read_count_threshold_;
    const size_t tail_threshold_;
};

}
}