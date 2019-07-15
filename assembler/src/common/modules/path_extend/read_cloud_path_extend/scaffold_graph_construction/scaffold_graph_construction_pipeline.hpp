#pragma once

#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/construction_callers.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/extender_searcher.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

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
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver::scaffold_graph_construction ScaffConParamsT;

    ScaffolderParams ConstructScaffolderParams(size_t length_threshold, const ScaffConParamsT &params,
                                               fragment_statistics::ClusterStatisticsExtractor primary_extractor) const;

    LongEdgePairGapCloserParams ConstructGapCloserParams(const ScaffConParamsT &params) const;

  private:
    ScaffolderParams::ScoreEstimationParams GetScoreEstimationParams(
        fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor,
        double score_percentile,
        double cluster_length_percentile,
        size_t block_length) const;

};

class ScaffoldGraphConstructionPipeline {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef shared_ptr<scaffold_graph::ScaffoldGraphConstructor> ConstructorPtr;

    ConstructorPtr initial_constructor_;
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> construction_stages_;
    vector<std::pair<shared_ptr<path_extend::scaffold_graph::ScaffoldGraph>, string>> intermediate_results_;
    const Graph &g_;
    const read_cloud::ScaffolderParams params_;

  public:
    ScaffoldGraphConstructionPipeline(ConstructorPtr initial_constructor, const Graph &g,
                                      const read_cloud::ScaffolderParams &params);

    void AddStage(shared_ptr<IterativeScaffoldGraphConstructorCaller> stage);

    void Run();

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraph> GetResult() const;

    vector<std::pair<shared_ptr<path_extend::scaffold_graph::ScaffoldGraph>, string>> GetIntermediateResults() const;
};

class ScaffoldGraphPipelineConstructor {
  protected:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> BarcodeIndexPtr;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor ScaffoldVertexExtractor;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigsT;

    const conj_graph_pack &gp_;
    const LibraryT lib_;
    const ReadCloudConfigsT read_cloud_configs_;

  public:
    ScaffoldGraphPipelineConstructor(const conj_graph_pack &gp, const LibraryT &lib, const ReadCloudConfigsT &configs);

    virtual ScaffoldGraphConstructionPipeline ConstructPipeline(const set<ScaffoldVertex> &scaffold_vertices) const = 0;

  protected:
    shared_ptr<ScaffoldVertexExtractor> ConstructSimpleEdgeIndex(const set<ScaffoldVertex> &scaffold_vertices,
                                                                 BarcodeIndexPtr barcode_extractor,
                                                                 const ScaffolderParams &params,
                                                                 size_t max_threads) const;
};

class BasicScaffoldGraphPipelineConstructor : public ScaffoldGraphPipelineConstructor {
  protected:
    using ScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using ScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;
    using ScaffoldGraphPipelineConstructor::gp_;
    using ScaffoldGraphPipelineConstructor::lib_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;
    size_t min_length_;

  public:
    BasicScaffoldGraphPipelineConstructor(const conj_graph_pack &gp,
                                          const LibraryT &lib,
                                          const ReadCloudConfigsT &configs,
                                          const ScaffoldingUniqueEdgeStorage &unique_storage,
                                          shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                          size_t max_threads_,
                                          size_t min_length_);

    ScaffoldGraphConstructionPipeline ConstructPipeline(const set<ScaffoldVertex> &scaffold_vertices) const override;
  private:
    virtual vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params,
        const set<ScaffoldVertex> &scaffold_vertices) const = 0;
};

class FullScaffoldGraphPipelineConstructor : public BasicScaffoldGraphPipelineConstructor {
    using BasicScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using BasicScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;
    using BasicScaffoldGraphPipelineConstructor::gp_;
    using BasicScaffoldGraphPipelineConstructor::lib_;
    using BasicScaffoldGraphPipelineConstructor::unique_storage_;
    using BasicScaffoldGraphPipelineConstructor::barcode_extractor_;
    using BasicScaffoldGraphPipelineConstructor::max_threads_;
    using BasicScaffoldGraphPipelineConstructor::min_length_;
    const ReadCloudSearchParameterPack search_parameter_pack_;

  public:
    FullScaffoldGraphPipelineConstructor(const conj_graph_pack &gp,
                                         const LibraryT &lib,
                                         const ReadCloudConfigsT &configs,
                                         const ScaffoldingUniqueEdgeStorage &unique_storage,
                                         shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                         size_t max_threads_,
                                         size_t min_length_,
                                         const ReadCloudSearchParameterPack &search_parameter_pack);

  private:
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params,
        const set<ScaffoldVertex> &scaffold_vertices) const override;
};

class GapScaffoldGraphPipelineConstructor : public ScaffoldGraphPipelineConstructor {
  protected:
    using ScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using ScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;
    using ScaffoldGraphPipelineConstructor::gp_;
    using ScaffoldGraphPipelineConstructor::lib_;
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;
    size_t min_length_;

  public:
    GapScaffoldGraphPipelineConstructor(const conj_graph_pack &gp,
                                        const LibraryT &lib,
                                        const ReadCloudConfigsT &configs,
                                        const ScaffoldingUniqueEdgeStorage &unique_storage,
                                        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                        size_t max_threads_,
                                        size_t min_length_);

    ScaffoldGraphConstructionPipeline ConstructPipeline(const set<ScaffoldVertex> &scaffold_vertices) const override;

  protected:
    shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetInitialConstructor(read_cloud::ScaffolderParams params,
                                                                               const set<ScaffoldVertex> &scaffold_vertices) const;

    virtual vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(
        read_cloud::ScaffolderParams params, const set<ScaffoldVertex> &scaffold_vertices) const = 0;
};

class MergingScaffoldGraphPipelineConstructor : public GapScaffoldGraphPipelineConstructor {
    using GapScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using GapScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;
    using GapScaffoldGraphPipelineConstructor::gp_;
    using GapScaffoldGraphPipelineConstructor::unique_storage_;
    using GapScaffoldGraphPipelineConstructor::barcode_extractor_;
    using GapScaffoldGraphPipelineConstructor::max_threads_;
    using GapScaffoldGraphPipelineConstructor::min_length_;

  public:
    MergingScaffoldGraphPipelineConstructor(const conj_graph_pack &gp,
                                            const LibraryT &lib,
                                            const ReadCloudConfigsT &configs,
                                            const ScaffoldingUniqueEdgeStorage &unique_storage,
                                            shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                            size_t max_threads_,
                                            size_t min_length_);

  protected:
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(ScaffolderParams params,
        const set<ScaffoldVertex> &scaffold_vertices) const override;
};

class BinningScaffoldGraphPipelineConstructor : public GapScaffoldGraphPipelineConstructor {
    using GapScaffoldGraphPipelineConstructor::ScaffoldGraph;
    using GapScaffoldGraphPipelineConstructor::ScaffoldVertex;
    using ScaffoldGraphPipelineConstructor::BarcodeIndexPtr;
    using GapScaffoldGraphPipelineConstructor::gp_;
    using GapScaffoldGraphPipelineConstructor::lib_;
    using GapScaffoldGraphPipelineConstructor::unique_storage_;
    using GapScaffoldGraphPipelineConstructor::barcode_extractor_;
    using GapScaffoldGraphPipelineConstructor::max_threads_;
    using GapScaffoldGraphPipelineConstructor::min_length_;

  public:
    BinningScaffoldGraphPipelineConstructor(const conj_graph_pack &gp,
                                            const LibraryT &lib,
                                            const ReadCloudConfigsT &configs,
                                            const ScaffoldingUniqueEdgeStorage &unique_storage,
                                            shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                            size_t max_threads,
                                            size_t min_length);

  private:
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(ScaffolderParams params,
        const set<ScaffoldVertex> &scaffold_vertices) const override;
};
}
}