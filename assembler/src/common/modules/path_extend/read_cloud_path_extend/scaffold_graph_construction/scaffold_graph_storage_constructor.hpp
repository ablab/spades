//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "scaffold_graph_construction_pipeline.hpp"
#include "extender_searcher.hpp"
#include "scaffold_graph_storage.hpp"
#include "modules/path_extend/pe_config_struct.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphStorageConstructor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef pe_config::ReadCloud ReadCloudConfigsT;

    ScaffoldGraphStorageConstructor(const ScaffoldingUniqueEdgeStorage &small_length_storage,
                                    const ScaffoldingUniqueEdgeStorage &large_length_storage,
                                    size_t small_length_threshold,
                                    size_t large_length_threshold,
                                    size_t max_threads,
                                    const LibraryT &lib,
                                    const ReadCloudConfigsT &configs,
                                    const ReadCloudSearchParameterPack &search_params,
                                    const std::string &scaffold_graph_path,
                                    const conj_graph_pack &gp);

    ScaffoldGraphStorage ConstructStorage() const;

    ScaffoldGraphStorage ConstructStorageFromPaths(const PathContainer &paths, bool scaffolding_mode) const;

  private:
    const ScaffoldingUniqueEdgeStorage &small_length_storage_;
    const ScaffoldingUniqueEdgeStorage &large_length_storage_;
    const size_t small_length_threshold_;
    const size_t large_length_threshold_;
    const size_t max_threads_;
    const LibraryT lib_;
    const ReadCloudConfigsT configs_;
    const ReadCloudSearchParameterPack search_params_;
    const std::string scaffold_graph_path_;
    const conj_graph_pack &gp_;
};

class ScaffoldGraphPolisherHelper {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef pe_config::ReadCloud CloudConfigT;
    ScaffoldGraphPolisherHelper(const Graph &g,
                                const debruijn_graph::Index &index,
                                const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                const CloudConfigT &cloud_configs,
                                size_t max_threads);

    ScaffoldGraph GetScaffoldGraphFromStorage(const ScaffoldGraphStorage &storage, bool path_scaffolding) const;

    void PrintScaffoldGraphReferenceInfo(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                         const ScaffoldingUniqueEdgeStorage &unique_storage,
                                         const std::string &path_to_reference) const;

  private:
    void GetGraphStorageReferenceInfo(const ScaffoldGraphStorage &storage, const std::string &path_to_reference) const;
    ScaffoldGraph ApplyRelativeThreshold(const ScaffoldGraph &graph,
                                         size_t unique_length_threshold,
                                         double relative_threshold) const;
    const Graph &g_;
    const debruijn_graph::Index &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    const CloudConfigT &cloud_configs_;
    const size_t max_threads_;
};

class CloudScaffoldGraphConstructor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::scaffolder::ScaffoldGraphConstructor ScaffoldGraphConstructor;
    typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef pe_config::ReadCloud ReadCloudConfigsT;
    typedef barcode_index::FrameBarcodeIndexInfoExtractor BarcodeExtractorT;

    CloudScaffoldGraphConstructor(size_t max_threads_,
                                  const debruijn_graph::conj_graph_pack &gp,
                                  const ScaffoldingUniqueEdgeStorage &unique_storage,
                                  const LibraryT &lib,
                                  const ReadCloudConfigsT &configs,
                                  const ReadCloudSearchParameterPack &search_parameter_pack,
                                  const std::string &debug_output_path,
                                  std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);
    ScaffoldGraph ConstructScaffoldGraphFromUniqueStorage(const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                          scaffold_graph_construction_pipeline_type::Type type) const;

    ScaffoldGraph ConstructScaffoldGraphFromPathContainer(const path_extend::PathContainer &paths,
                                                          size_t min_length, bool scaffolding_mode) const;

    ScaffoldGraph ConstructScaffoldGraphFromVertices(const std::set<ScaffoldVertex> &scaffold_vertices,
                                                     const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                     size_t min_length,
                                                     scaffold_graph_construction_pipeline_type::Type type) const;

  private:
    const size_t max_threads_;
    const debruijn_graph::conj_graph_pack &gp_;
    const ScaffoldingUniqueEdgeStorage &unique_storage_;
    const LibraryT lib_;
    const ReadCloudConfigsT &configs_;
    const ReadCloudSearchParameterPack search_parameter_pack_;
    const std::string debug_output_path_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

    DECL_LOGGER("CloudScaffoldGraphConstructor");
};
}
}