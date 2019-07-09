#pragma once

#include "scaffold_graph_construction_pipeline.hpp"
#include "scaffold_graph_storage.hpp"

namespace path_extend {
class ScaffoldGraphStorageConstructor {
 public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
 private:
    const size_t small_length_threshold_;
    const size_t large_length_threshold_;
    const LibraryT lib_;
    const conj_graph_pack &gp_;

 public:
    ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                    size_t large_length_threshold_,
                                    const LibraryT &lib,
                                    const conj_graph_pack &gp_);

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
                                         const std::string &path_to_reference,
                                         size_t length_threshold) const;

 private:
    void GetGraphStorageReferenceInfo(const ScaffoldGraphStorage &storage, const std::string &path_to_reference) const;

    ScaffoldGraph ApplyRelativeThreshold(const ScaffoldGraph &graph, size_t unique_length_threshold, double relative_threshold) const;
};

class CloudScaffoldGraphConstructor {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::scaffold_graph::ScaffoldGraphConstructor ScaffoldGraphConstructor;
    typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
    typedef io::SequencingLibrary<debruijn_graph::config::LibraryData> LibraryT;
    typedef barcode_index::FrameBarcodeIndexInfoExtractor BarcodeExtractorT;
 private:
    const size_t max_threads_;
    const debruijn_graph::conj_graph_pack &gp_;
    const LibraryT lib_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;


 public:
    CloudScaffoldGraphConstructor(size_t max_threads_,
                                  const debruijn_graph::conj_graph_pack &gp,
                                  const LibraryT &lib,
                                  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);
    ScaffoldGraph ConstructScaffoldGraphFromMinLength(size_t min_length,
                                                      scaffold_graph_construction_pipeline_type::Type type) const;

    ScaffoldGraph ConstructScaffoldGraphFromPathContainer(const path_extend::PathContainer &paths,
                                                          size_t min_length, bool scaffolding_mode) const;

    ScaffoldGraph ConstructScaffoldGraphFromVertices(const set<ScaffoldVertex> &scaffold_vertices,
                                                     const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                     size_t min_length,
                                                     scaffold_graph_construction_pipeline_type::Type type) const;

    ScaffoldingUniqueEdgeStorage ConstructUniqueStorage(size_t min_length) const;

    DECL_LOGGER("CloudScaffoldGraphConstructor");
};
}