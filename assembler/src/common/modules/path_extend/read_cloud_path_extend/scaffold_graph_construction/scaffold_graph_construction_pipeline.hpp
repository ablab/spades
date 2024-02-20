//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "construction_callers.hpp"
#include "extender_searcher.hpp"
#include "read_cloud_connection_conditions.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_vertex.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/cluster_distribution_estimator.hpp"

namespace path_extend {
namespace read_cloud {

namespace scaffold_graph_construction_pipeline_type {
enum Type {
  Basic,
  Scaffolding,
  Binning,
};
}

class ScaffolderParamsConstructor {
  public:
    typedef pe_config::ReadCloud::scaffold_graph_construction ScaffConParamsT;

    ScaffolderParams ConstructScaffolderParams(size_t length_threshold, const ScaffConParamsT &params,
                                               fragment_statistics::ClusterStatisticsExtractor primary_extractor) const;

    LongEdgePairGapCloserParams ConstructGapCloserParams(const ScaffConParamsT &params) const;

  private:
    ScaffolderParams::ScoreEstimationParams GetScoreEstimationParams(
        fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor,
        double score_percentile,
        double cluster_length_percentile,
        size_t block_length, size_t initial_distance) const;

};

class ScaffoldGraphConstructionPipeline {
  public:
    typedef graph_pack::GraphPack GraphPack;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef std::shared_ptr<scaffolder::ScaffoldGraphConstructor> ConstructorPtr;
    typedef std::pair<std::shared_ptr<scaffold_graph::ScaffoldGraph>, string> ResultT;

    ScaffoldGraphConstructionPipeline(ConstructorPtr initial_constructor,
                                      const read_cloud::ScaffolderParams &params,
                                      const Graph &g);

    void AddStage(std::shared_ptr<IterativeScaffoldGraphConstructorCaller> stage);
    void Run();
    std::shared_ptr<scaffold_graph::ScaffoldGraph> GetResult() const;

    std::vector<ResultT> GetIntermediateResults() const;

  private:
    ConstructorPtr initial_constructor_;
    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> construction_stages_;
    std::vector<ResultT> intermediate_results_;
    const read_cloud::ScaffolderParams params_;
    const Graph &g_;
};

class ScaffoldGraphPipelineConstructor {
  public:
    typedef graph_pack::GraphPack GraphPack;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> BarcodeIndexPtr;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor ScaffoldVertexExtractor;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef pe_config::ReadCloud ReadCloudConfigsT;

    ScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs, const Graph &g);
    virtual ~ScaffoldGraphPipelineConstructor() = default;

    virtual ScaffoldGraphConstructionPipeline ConstructPipeline(const std::set<ScaffoldVertex> &scaffold_vertices) const = 0;

  protected:
    std::shared_ptr<ScaffoldVertexExtractor> ConstructSimpleEdgeIndex(const std::set<ScaffoldVertex> &scaffold_vertices,
                                                                      BarcodeIndexPtr barcode_extractor,
                                                                      const ScaffolderParams &params,
                                                                      size_t max_threads) const;
    const ReadCloudConfigsT read_cloud_configs_;
    const Graph &g_;
};

class BasicScaffoldGraphPipelineConstructor : public ScaffoldGraphPipelineConstructor {
  public:
    BasicScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs,
                                          const graph_pack::GraphPack &gp,
                                          const LibraryT &lib,
                                          const ScaffoldingUniqueEdgeStorage &unique_storage,
                                          std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                          size_t max_threads,
                                          size_t min_length);

    virtual ~BasicScaffoldGraphPipelineConstructor() = default;

    ScaffoldGraphConstructionPipeline ConstructPipeline(const std::set<ScaffoldVertex> &scaffold_vertices) const override;

  protected:
    using ScaffoldGraphPipelineConstructor::GraphPack;
    using ScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using ScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;

    const graph_pack::GraphPack &gp_;
    const LibraryT lib_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;
    size_t min_length_;

  private:
    virtual std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params,
        const std::set<ScaffoldVertex> &scaffold_vertices) const = 0;
};

class FullScaffoldGraphPipelineConstructor : public BasicScaffoldGraphPipelineConstructor {
  public:
    using BasicScaffoldGraphPipelineConstructor::GraphPack;
    using BasicScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using BasicScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;

    FullScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs,
                                         const graph_pack::GraphPack &gp,
                                         const LibraryT &lib,
                                         const ScaffoldingUniqueEdgeStorage &unique_storage,
                                         std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                         size_t max_threads,
                                         size_t min_length,
                                         const ReadCloudSearchParameterPack &search_parameter_pack);

  private:
    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params,
        const std::set<ScaffoldVertex> &scaffold_vertices) const override;

    using BasicScaffoldGraphPipelineConstructor::gp_;
    using BasicScaffoldGraphPipelineConstructor::lib_;
    using BasicScaffoldGraphPipelineConstructor::unique_storage_;
    using BasicScaffoldGraphPipelineConstructor::barcode_extractor_;
    using BasicScaffoldGraphPipelineConstructor::max_threads_;
    using BasicScaffoldGraphPipelineConstructor::min_length_;
    const ReadCloudSearchParameterPack search_parameter_pack_;
};

class GapScaffoldGraphPipelineConstructor : public ScaffoldGraphPipelineConstructor {
  public:
    using ScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using ScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;

    GapScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs,
                                        const Graph &g,
                                        const LibraryT &lib,
                                        const ScaffoldingUniqueEdgeStorage &unique_storage,
                                        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                        size_t max_threads,
                                        size_t min_length);

    virtual ~GapScaffoldGraphPipelineConstructor() = default;

    ScaffoldGraphConstructionPipeline ConstructPipeline(const std::set<ScaffoldVertex> &scaffold_vertices) const override;

  protected:
    std::shared_ptr<scaffolder::ScaffoldGraphConstructor> GetInitialConstructor(
        read_cloud::ScaffolderParams params,
        const std::set<ScaffoldVertex> &scaffold_vertices) const;

    virtual std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params, const std::set<ScaffoldVertex> &scaffold_vertices) const = 0;

    using ScaffoldGraphPipelineConstructor::g_;
    const LibraryT &lib_;
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;
    size_t min_length_;
};

class MergingScaffoldGraphPipelineConstructor : public GapScaffoldGraphPipelineConstructor {
  public:
    using GapScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using GapScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;

    MergingScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs,
                                            const Graph &g,
                                            const LibraryT &lib,
                                            const ScaffoldingUniqueEdgeStorage &unique_storage,
                                            std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                            size_t max_threads,
                                            size_t min_length);

  protected:
    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        ScaffolderParams params,
        const std::set<ScaffoldVertex> &scaffold_vertices) const override;

    using GapScaffoldGraphPipelineConstructor::g_;
    using GapScaffoldGraphPipelineConstructor::unique_storage_;
    using GapScaffoldGraphPipelineConstructor::barcode_extractor_;
    using GapScaffoldGraphPipelineConstructor::max_threads_;
    using GapScaffoldGraphPipelineConstructor::min_length_;
};

class BinningScaffoldGraphPipelineConstructor : public GapScaffoldGraphPipelineConstructor {
  public:
    using GapScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using GapScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;

    BinningScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs,
                                            const Graph &g,
                                            const LibraryT &lib,
                                            const ScaffoldingUniqueEdgeStorage &unique_storage,
                                            std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                            size_t max_threads,
                                            size_t min_length);

  protected:
    using GapScaffoldGraphPipelineConstructor::g_;
    using GapScaffoldGraphPipelineConstructor::lib_;
    using GapScaffoldGraphPipelineConstructor::unique_storage_;
    using GapScaffoldGraphPipelineConstructor::barcode_extractor_;
    using GapScaffoldGraphPipelineConstructor::max_threads_;
    using GapScaffoldGraphPipelineConstructor::min_length_;

  private:
    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        ScaffolderParams params,
        const std::set<ScaffoldVertex> &scaffold_vertices) const override;
};
}
}