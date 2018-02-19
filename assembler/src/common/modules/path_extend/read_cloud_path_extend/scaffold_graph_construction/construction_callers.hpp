#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"

namespace path_extend {

struct ScaffolderParams {
  size_t length_threshold_;
  size_t tail_threshold_;
  size_t count_threshold_;
  double vertex_multiplier_;
  double connection_score_threshold_;
  double relative_coverage_threshold_;
  size_t connection_length_threshold_;
  size_t connection_count_threshold_;
  size_t initial_distance_;
  double split_procedure_strictness_;
  size_t transitive_distance_threshold_;
  size_t min_length_for_barcode_collection_;

  ScaffolderParams(size_t length_threshold_,
                   size_t tail_threshold_,
                   size_t count_threshold_,
                   double vertex_multiplier_,
                   double connection_barcode_threshold_,
                   double relative_coverage_threshold_,
                   size_t connection_length_threshold_,
                   size_t connection_count_threshold_,
                   size_t initial_distance_,
                   double split_procedure_strictness_,
                   size_t transitive_distance_threshold_,
                   size_t min_length_for_barcode_collection);
};

class IterativeScaffoldGraphConstructorCaller {
 protected:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 public:
    virtual shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                         const ScaffoldGraph& scaffold_graph) const = 0;

    virtual ~IterativeScaffoldGraphConstructorCaller() = default;
};
class BarcodeScoreConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;

 public:
    BarcodeScoreConstructorCaller(const Graph& g_,
                                  shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                  size_t max_threads_);
    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                 const ScaffoldGraph& scaffold_graph) const override;
};
class BarcodeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    size_t max_threads;

 public:
    BarcodeConnectionConstructorCaller(const Graph& g_,
                                       shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                       shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                       const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
                                       size_t max_threads);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                 const ScaffoldGraph& scaffold_graph) const override;
};
class CompositeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const debruijn_graph::conj_graph_pack& gp_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const size_t max_threads_;
    const bool scaffolding_mode_;

 public:
    CompositeConnectionConstructorCaller(const debruijn_graph::conj_graph_pack &gp_,
                                         shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                         shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                         const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_,
                                         const size_t max_threads_, bool scaffolding_mode);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                 const ScaffoldGraph& scaffold_graph) const override;
};
class EdgeSplitConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    std::size_t max_threads_;

 public:
    EdgeSplitConstructorCaller(const Graph& g_,
                               shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                               std::size_t max_threads_);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                 const ScaffoldGraph& scaffold_graph) const override;
};
class TransitiveConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    std::size_t max_threads_;

 public:
    TransitiveConstructorCaller(const Graph& g_,
                                std::size_t max_threads_);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const path_extend::ScaffolderParams& params,
                                                                                                 const ScaffoldGraph& scaffold_graph) const override;
};
} //path_extend