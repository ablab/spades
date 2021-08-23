//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "extender_searcher.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "modules/path_extend/path_extender.hpp"
#include "modules/path_extend/pipeline/launch_support.hpp"
#include "modules/path_extend/extension_chooser.hpp"
#include "modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/pair_entry_processors.hpp"

namespace path_extend {
namespace read_cloud {

//Same as AssemblyGraphConnectionCondition, but stops after reaching unique edges.
class AssemblyGraphUniqueConnectionCondition : public AssemblyGraphConnectionCondition {
  public:
    AssemblyGraphUniqueConnectionCondition(const Graph &g,
                                           size_t max_connection_length,
                                           const ScaffoldingUniqueEdgeStorage &unique_edges);
    std::map<EdgeId, double> ConnectedWith(EdgeId e) const override;
    bool IsLast() const override;

  private:
    using AssemblyGraphConnectionCondition::g_;
    using AssemblyGraphConnectionCondition::interesting_edge_set_;
    using AssemblyGraphConnectionCondition::max_connection_length_;
    //fixme duplication with interesting edges, needed to pass to dijkstra
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
};

class ScaffoldEdgePredicate : public func::AbstractPredicate<const scaffold_graph::ScaffoldGraph::ScaffoldEdge &> {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

    virtual ~ScaffoldEdgePredicate() = default;
};

struct ReadCloudMiddleDijkstraParams {
  ReadCloudMiddleDijkstraParams(size_t count_threshold_,
                                size_t tail_threshold_,
                                size_t distance_,
                                const LongEdgePairGapCloserParams &edge_pair_gap_closer_params_);

  const size_t count_threshold_;
  const size_t tail_threshold_;
  const size_t distance_;
  const LongEdgePairGapCloserParams edge_pair_gap_closer_params_;
};

class ReadCloudMiddleDijkstraPredicate : public ScaffoldEdgePredicate {
  public:
    using ScaffoldEdgePredicate::ScaffoldEdge;
    ReadCloudMiddleDijkstraPredicate(const Graph &g,
                                     const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                     std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
                                     std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                     const ReadCloudMiddleDijkstraParams &params);
    bool Check(const ScaffoldEdge &scaffold_edge) const override;

  private:
    const Graph &g;
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
    const std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor_;
    const std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    const ReadCloudMiddleDijkstraParams params_;
    DECL_LOGGER("ReadCloudMiddleDijkstraPredicate");
};

//class CompositeConnectionPredicate : public ScaffoldEdgePredicate {
//  public:
//    using ScaffoldEdgePredicate::ScaffoldEdge;
//    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
//
//    CompositeConnectionPredicate(const conj_graph_pack &gp,
//                                 std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
//                                 std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
//                                 const ScaffoldingUniqueEdgeStorage &unique_storage,
//                                 size_t length_bound,
//                                 const ReadCloudSearchParameterPack &search_parameter_pack,
//                                 const LongEdgePairGapCloserParams &predicate_params,
//                                 bool scaffolding_mode);
//
//    bool Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &scaffold_edge) const override;
//
//  private:
//    std::shared_ptr<scaffolder::ScaffoldVertexPredicate> ConstructScaffoldVertexPredicate(
//        const ScaffoldVertex &start, const ScaffoldVertex &end,
//        std::shared_ptr<PairEntryProcessor> entry_processor) const;
//
//    const conj_graph_pack &gp_;
//    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor_;
//    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
//    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
//    const size_t length_bound_;
//    const ReadCloudSearchParameterPack search_parameter_pack_;
//    const LongEdgePairGapCloserParams predicate_params_;
//    bool scaffolding_mode_;
//
//    DECL_LOGGER("CompositeConnectionPredicate");
//};

class EdgeSplitPredicate : public ScaffoldEdgePredicate {
  public:
    using ScaffoldEdgePredicate::ScaffoldEdge;
    typedef barcode_index::BarcodeId BarcodeId;
    EdgeSplitPredicate(const Graph &g_,
                       std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
                       size_t count_threshold_,
                       double strictness);

    bool Check(const ScaffoldEdge &scaffold_edge) const override;

  private:
    bool CheckOrderingForThreeSegments(const barcode_index::SimpleVertexEntry &first,
                                       const barcode_index::SimpleVertexEntry &second,
                                       const barcode_index::SimpleVertexEntry &third, double strictness) const;

    bool CheckOrderingForFourSegments(const barcode_index::SimpleVertexEntry &first,
                                      const barcode_index::SimpleVertexEntry &second,
                                      const barcode_index::SimpleVertexEntry &third,
                                      const barcode_index::SimpleVertexEntry &fourth) const;

    const Graph &g_;
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    const size_t count_threshold_;
    const double strictness_;

    DECL_LOGGER("EdgeSplitPredicate");
};

class EdgeInTheMiddlePredicate {
  public:
    typedef barcode_index::BarcodeId BarcodeId;

    EdgeInTheMiddlePredicate(const Graph &g_,
                             const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                             size_t count_threshold,
                             double shared_fraction_threshold);

    bool IsCorrectOrdering(const EdgeId &first, const EdgeId &second, const EdgeId &third);

  private:
    const Graph &g_;
    const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_;
    const size_t count_threshold_;
    const double shared_fraction_threshold_;
    DECL_LOGGER("EdgeInTheMiddlePredicate");
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