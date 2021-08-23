//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "extender_searcher.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"

#include <memory>

namespace path_extend {
namespace read_cloud {

struct ScaffolderParams {
  struct ScoreEstimationParams {
    double score_percentile_;
    size_t max_cluster_gap_;
    size_t training_edge_length_threshold_;

    ScoreEstimationParams(double score_percentile, size_t max_distance, size_t training_edge_length_threshold);
  };

  ScaffolderParams(size_t length_threshold,
                   size_t tail_threshold,
                   size_t count_threshold,
                   size_t connection_length_threshold,
                   size_t connection_count_threshold,
                   size_t initial_distance,
                   double split_procedure_strictness,
                   size_t transitive_distance_threshold,
                   size_t min_length_for_barcode_collection,
                   const LongEdgePairGapCloserParams &gap_closer_params,
                   const ScoreEstimationParams &score_estimation_params);

  size_t length_threshold_;
  size_t tail_threshold_;
  size_t count_threshold_;
  size_t connection_length_threshold_;
  size_t connection_count_threshold_;
  size_t initial_distance_;
  double split_procedure_strictness_;
  size_t transitive_distance_threshold_;
  size_t min_length_for_barcode_collection_;
  LongEdgePairGapCloserParams gap_closer_params_;
  ScoreEstimationParams score_estimation_params_;
};

/** Interface for scaffold graph factory. Implementations use initial scaffold graph to construct the new one.
 */
class IterativeScaffoldGraphConstructorCaller {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef std::shared_ptr<scaffolder::ScaffoldGraphConstructor> GraphConstructor;
    explicit IterativeScaffoldGraphConstructorCaller(const string &name);
    /**
    * @param scaffold_graph Original scaffold graph which is used to construct the new one.
     * Most constructors remove subset of edges of the original graph based on some predicate.
    */
    virtual GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const = 0;

    virtual ~IterativeScaffoldGraphConstructorCaller() = default;
    string getName() const;

  private:
    string name_;
};
class ScoreFunctionConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
  public:
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
    typedef std::shared_ptr<path_extend::read_cloud::ScaffoldEdgeScoreFunction> ScoreFunction;
    ScoreFunctionConstructorCaller(const string &name, const Graph &g, size_t max_threads);
    GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const override;
  protected:
    const Graph &g_;
    size_t max_threads_;
  private:
    virtual ScoreFunction ConstructScoreFunction() const = 0;
    virtual double ConstructScoreThreshold() const = 0;
};
/** ConstructorCaller that removes edges from scaffold graph based on intersection between
 * sets of barcodes from two long edges
 */
class BarcodeScoreConstructorCaller : public ScoreFunctionConstructorCaller {
  public:
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
    using ScoreFunctionConstructorCaller::ScoreFunction;
    BarcodeScoreConstructorCaller(const Graph &g,
                                  size_t max_threads,
                                  std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor,
                                  std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor,
                                  const ScaffolderParams &params);
  private:
    ScoreFunction ConstructScoreFunction() const override;
    double ConstructScoreThreshold() const override;

    using ScoreFunctionConstructorCaller::g_;
    using ScoreFunctionConstructorCaller::max_threads_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor_;
    std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    ScaffolderParams params_;
};
/** ConstructorCaller that ascertains existence of barcode-supported short-edge path between two long edges
 */
class BarcodeConnectionConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
  public:
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
    BarcodeConnectionConstructorCaller(const Graph &g,
                                       std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                       std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                       const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
                                       const ScaffolderParams &params,
                                       size_t max_threads);

    GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const override;
  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
    ScaffolderParams params_;
    size_t max_threads_;
};
/** ConstructorCaller that ascertains existence of short-edge path between two long edges that is
 * supported by barcodes and paired info
 */
//class CompositeConnectionConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
//  public:
//    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
//    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
//    typedef std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> MainBarcodeIndexPtr;
//    typedef std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> ScaffoldBarcodeIndexPtr;
//    typedef pe_config::ReadCloud::scaffold_graph_construction ScaffConConfigs;
//
//    CompositeConnectionConstructorCaller(const debruijn_graph::conj_graph_pack &gp,
//                                         MainBarcodeIndexPtr main_extractor,
//                                         ScaffoldBarcodeIndexPtr barcode_extractor,
//                                         const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
//                                         const ReadCloudSearchParameterPack &search_parameter_pack,
//                                         const ScaffConConfigs &scaff_con_configs,
//                                         const ScaffolderParams &params,
//                                         size_t max_threads,
//                                         bool scaffolding_mode);
//
//    GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const override;
//
//  private:
//    const debruijn_graph::conj_graph_pack &gp_;
//    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
//    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
//    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
//    const ReadCloudSearchParameterPack search_parameter_pack_;
//    const ScaffConConfigs &scaff_con_configs_;
//    ScaffolderParams params_;
//    size_t max_threads_;
//    bool scaffolding_mode_;
//
//    DECL_LOGGER("CompositeConnectionConstructorCaller");
//};
/** ConstructorCaller that filters conjugate transitions by selecting the best orientation of two edges in the transition.
 */
class EdgeSplitConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
  public:
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
    EdgeSplitConstructorCaller(const Graph &g,
                               std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
                               const ScaffolderParams &params,
                               size_t max_threads);

    GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const override;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    ScaffolderParams params_;
    size_t max_threads_;
};
/** ConstructorCaller that filters transitive connections.
 */
class TransitiveConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
  public:
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    using IterativeScaffoldGraphConstructorCaller::GraphConstructor;
    TransitiveConstructorCaller(const Graph &g,
                                size_t max_threads,
                                size_t transitive_distance_threshold);

    GraphConstructor GetScaffoldGraphConstuctor(const ScaffoldGraph &scaffold_graph) const override;

  private:
    const Graph &g_;
    size_t max_threads_;
    size_t transitive_distance_threshold_;
};

} //path_extend
} //read_cloud