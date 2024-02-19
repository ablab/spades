#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"

namespace path_extend {

struct ScaffolderParams {
  struct ScoreEstimationParams {
    double score_percentile_;
    size_t max_cluster_gap_;
    size_t training_edge_length_threshold_;

    ScoreEstimationParams(double score_percentile_, size_t max_distance_, size_t training_edge_length_threshold_);
  };

  size_t length_threshold_;
  size_t tail_threshold_;
  size_t count_threshold_;
  size_t connection_length_threshold_;
  size_t connection_count_threshold_;
  size_t initial_distance_;
  double split_procedure_strictness_;
  size_t transitive_distance_threshold_;
  size_t min_length_for_barcode_collection_;
  ScoreEstimationParams score_estimation_params_;

  ScaffolderParams(size_t length_threshold_,
                   size_t tail_threshold_,
                   size_t count_threshold_,
                   size_t connection_length_threshold_,
                   size_t connection_count_threshold_,
                   size_t initial_distance_,
                   double split_procedure_strictness_,
                   size_t transitive_distance_threshold_,
                   size_t min_length_for_barcode_collection,
                   const ScoreEstimationParams& score_estimation_params);
};

/** Interface for scaffold graph factory. Implementations use initial scaffold graph to construct new one.
 */
class IterativeScaffoldGraphConstructorCaller {
 protected:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

 private:
    string name_;
 public:
    IterativeScaffoldGraphConstructorCaller(const string &name);
    /**
    * @param params Scaffold graph construction configs
    * @param scaffold_graph Original scaffold graph which is used to construct new one.
     * Most constructors remove subset of edges of original graph based on a predicate.
    */
    virtual shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const = 0;

    virtual ~IterativeScaffoldGraphConstructorCaller() = default;

    string getName() const;
};
/** ConstructorCaller that removes edges from scaffold graph based on containment index between
 * sets of barcodes of two long edges
 */
class BarcodeScoreConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor_;
    shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;

 public:
    BarcodeScoreConstructorCaller(const Graph& g_,
                                  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor,
                                  shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                  size_t max_threads_);
    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const override;
};
/** ConstructorCaller that ascertains existence of cloud-supported path between two long edges
 */
class BarcodeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    size_t max_threads_;

 public:
    BarcodeConnectionConstructorCaller(const Graph& g_,
                                       shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                       shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                       const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
                                       size_t max_threads);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const override;
};
/** ConstructorCaller that ascertains existence of path between two edges that is supported by clouds and paired info
 */
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

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const override;
};
/** ConstructorCaller that filters conjugate connections by selecting best orientation of two edges in connection.
 */
class EdgeSplitConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_;
    std::size_t max_threads_;

 public:
    EdgeSplitConstructorCaller(const Graph& g_,
                               shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                               std::size_t max_threads_);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const override;
};
/** ConstructorCaller that filters transitive connections.
 */
class TransitiveConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
    using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
    const Graph& g_;
    std::size_t max_threads_;

 public:
    TransitiveConstructorCaller(const Graph& g_,
                                std::size_t max_threads_);

    shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(
        const path_extend::ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const override;
};
} //path_extend