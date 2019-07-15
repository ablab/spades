#pragma once

#include "scaffold_graph_construction_pipeline.hpp"
#include "extender_searcher.hpp"
#include "scaffold_graph_storage.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphStorageConstructor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigsT;
  private:
    const ScaffoldingUniqueEdgeStorage &small_length_storage_;
    const ScaffoldingUniqueEdgeStorage &large_length_storage_;
    const size_t small_length_threshold_;
    const size_t large_length_threshold_;
    const size_t max_threads_;
    const LibraryT lib_;
    const ReadCloudConfigsT configs_;
    const ReadCloudSearchParameterPack search_params_;
    const conj_graph_pack &gp_;

  public:
    ScaffoldGraphStorageConstructor(const ScaffoldingUniqueEdgeStorage &small_length_storage,
                                    const ScaffoldingUniqueEdgeStorage &large_length_storage,
                                    size_t small_length_threshold,
                                    size_t large_length_threshold,
                                    size_t max_threads,
                                    const LibraryT &lib,
                                    const ReadCloudConfigsT &configs,
                                    const ReadCloudSearchParameterPack &search_params,
                                    const conj_graph_pack &gp);

    ScaffoldGraphStorage ConstructStorage() const;

    ScaffoldGraphStorage ConstructStorageFromPaths(const PathContainer &paths, bool scaffolding_mode) const;
};

class ScaffoldGraphPolisherHelper {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef config::debruijn_config::read_cloud_resolver CloudConfigT;
  private:
    const conj_graph_pack &gp_;
    const CloudConfigT &cloud_configs_;
    const size_t max_threads_;
  public:
    ScaffoldGraphPolisherHelper(const conj_graph_pack &gp, const CloudConfigT &cloud_configs, size_t max_threads);

    ScaffoldGraph GetScaffoldGraphFromStorage(const ScaffoldGraphStorage &storage, bool path_scaffolding) const;

    void PrintScaffoldGraphReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                         const std::string &path_to_reference) const;

  private:
    void GetGraphStorageReferenceInfo(const ScaffoldGraphStorage &storage, const std::string &path_to_reference) const;

    ScaffoldGraph ApplyRelativeThreshold(const ScaffoldGraph &graph,
                                         size_t unique_length_threshold,
                                         double relative_threshold) const;
};

class CloudScaffoldGraphConstructor {
  public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::scaffold_graph::ScaffoldGraphConstructor ScaffoldGraphConstructor;
    typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigsT;
    typedef barcode_index::FrameBarcodeIndexInfoExtractor BarcodeExtractorT;
  private:
    const size_t max_threads_;
    const debruijn_graph::conj_graph_pack &gp_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
    const LibraryT lib_;
    const ReadCloudConfigsT &configs_;
    const ReadCloudSearchParameterPack search_parameter_pack_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

  public:
    CloudScaffoldGraphConstructor(size_t max_threads_,
                                  const debruijn_graph::conj_graph_pack &gp,
                                  const ScaffoldingUniqueEdgeStorage &unique_storage,
                                  const LibraryT &lib,
                                  const ReadCloudConfigsT &configs,
                                  const ReadCloudSearchParameterPack search_parameter_pack,
                                  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);
    ScaffoldGraph ConstructScaffoldGraphFromMinLength(const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                      scaffold_graph_construction_pipeline_type::Type type) const;

    ScaffoldGraph ConstructScaffoldGraphFromPathContainer(const path_extend::PathContainer &paths,
                                                          size_t min_length, bool scaffolding_mode) const;

    ScaffoldGraph ConstructScaffoldGraphFromVertices(const set<ScaffoldVertex> &scaffold_vertices,
                                                     const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                     size_t min_length,
                                                     scaffold_graph_construction_pipeline_type::Type type) const;

    DECL_LOGGER("CloudScaffoldGraphConstructor");
};
}
}