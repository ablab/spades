#pragma once

#include "scaffold_graph_construction_pipeline.hpp"
#include "scaffold_graph_storage.hpp"

namespace path_extend {
class ScaffoldGraphStorageConstructor {
 public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const size_t small_length_threshold_;
    const size_t large_length_threshold_;
    const conj_graph_pack &gp_;

 public:
    ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                    size_t large_length_threshold_,
                                    const conj_graph_pack &gp_);

    ScaffoldGraphStorage ConstructStorageFromGraph() const;

    ScaffoldGraphStorage ConstructStorageFromPaths(const PathContainer &paths, bool scaffolding_mode) const;
};

class ScaffoldGraphPolisherLauncher {
 public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const conj_graph_pack &gp_;
 public:
    ScaffoldGraphPolisherLauncher(const conj_graph_pack &gp_);

    ScaffoldGraph GetScaffoldGraphFromStorage(const ScaffoldGraphStorage &storage, bool path_scaffolding) const;

 private:
    void GetGraphStorageReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &small_scaffold_graph,
                                      const path_extend::scaffold_graph::ScaffoldGraph &large_scaffold_graph,
                                      const debruijn_graph::conj_graph_pack &graph_pack) const;

    void PrintScaffoldGraphReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                         const debruijn_graph::conj_graph_pack &graph_pack,
                                         size_t length_threshold) const;
};

class CloudScaffoldGraphConstructor {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::scaffold_graph::ScaffoldGraphConstructor ScaffoldGraphConstructor;
    typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
 private:
    const size_t max_threads_;
    const debruijn_graph::conj_graph_pack &gp_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

 public:
    CloudScaffoldGraphConstructor(size_t max_threads_,
                                  const debruijn_graph::conj_graph_pack &gp,
                                  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);
    ScaffoldGraph ConstructScaffoldGraphFromMinLength(size_t min_length,
                                                      scaffold_graph_construction_pipeline_type::Type type) const;

    ScaffoldGraph ConstructScaffoldGraphFromPathContainer(const path_extend::PathContainer &paths,
                                                          size_t min_length, bool scaffolding_mode) const;

    ScaffoldGraph ConstructScaffoldGraphFromVertices(const set<ScaffoldVertex> &scaffold_vertices,
                                                     const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                     size_t min_length,
                                                     scaffold_graph_construction_pipeline_type::Type type) const;

 private:
    ScaffoldingUniqueEdgeStorage ConstructUniqueStorage(size_t min_length) const;
};
}