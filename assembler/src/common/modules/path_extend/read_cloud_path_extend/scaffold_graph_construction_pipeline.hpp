#pragma once

#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "read_cloud_connection_conditions.hpp"
#include "scaffold_graph_storage.hpp"

namespace path_extend {
    struct ScaffolderParams {
      const size_t length_threshold_;
      const size_t tail_threshold_;
      const size_t count_threshold_;
      const double score_threshold_;
      const double connection_score_threshold_;
      const size_t connection_length_threshold_;
      const size_t connection_count_threshold_;
      const size_t initial_distance_;
      const double split_procedure_strictness_;
      const size_t transitive_distance_threshold_;

      ScaffolderParams(size_t length_threshold_,
                       size_t tail_threshold_,
                       size_t count_threshold_,
                       double score_threshold_,
                       double connection_barcode_threshold_,
                       size_t connection_length_threshold_,
                       size_t connection_count_threshold_,
                       size_t initial_distance_,
                       double split_procedure_strictness_,
                       size_t transitive_distance_threshold_);
    };

    class ScaffolderParamsConstructor {
     public:
        ScaffolderParams ConstructScaffolderParams(size_t length_threshold);
    };

    class IterativeScaffoldGraphConstructorCaller {
     protected:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
     public:
        virtual shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                                const ScaffoldGraph& scaffold_graph) const = 0;

        virtual ~IterativeScaffoldGraphConstructorCaller() = default;
    };

    class BarcodeScoreConstructorCaller : public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        size_t max_threads_;

     public:
        BarcodeScoreConstructorCaller(const Graph& g_,
                                      const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                      size_t max_threads_);
        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class BarcodeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
        size_t max_threads;

     public:
        BarcodeConnectionConstructorCaller(const Graph& g_,
                                           const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                           const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                           size_t max_threads);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class CompositeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const conj_graph_pack& gp_;
        const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
        const size_t max_threads_;

     public:
        CompositeConnectionConstructorCaller(const conj_graph_pack &gp_,
                                             const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                             const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                             const size_t max_threads_);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class EdgeSplitConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        size_t max_threads_;

     public:
        EdgeSplitConstructorCaller(const Graph& g_,
                                   const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                   size_t max_threads_);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class TransitiveConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
        size_t max_threads_;

     public:
        TransitiveConstructorCaller(const Graph& g_,
                                    const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                    size_t max_threads_);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class CloudScaffoldGraphConstuctor {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef scaffold_graph::ScaffoldGraphConstructor ScaffoldGraphConstructor;
        typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
     private:
        const size_t max_threads_;
        const size_t full_pipeline_length_;
        const conj_graph_pack& gp_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;

     public:
        CloudScaffoldGraphConstuctor(size_t max_threads_,
                                     size_t full_pipeline_length_threshold,
                                     const conj_graph_pack& gp,
                                     const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor);
        ScaffoldGraph ConstructScaffoldGraphFromLength(size_t min_length) const;

        ScaffoldGraph ConstructScaffoldGraphFromStorage(const ScaffolderParams& params, const ScaffoldingUniqueEdgeStorage& unique_storage) const;

     private:

        vector<shared_ptr<ConnectionCondition>> GetGraphConnectionConditions(const ScaffolderParams& params,
                                                                             const ScaffoldingUniqueEdgeStorage& unique_storage) const;
    };


    class ScaffoldGraphStorageConstructor {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
     private:
        const size_t small_length_threshold_;
        const size_t large_length_threshold_;
        const conj_graph_pack& gp_;

     public:
        ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                        size_t large_length_threshold_,
                                        const conj_graph_pack& gp_);

        ScaffoldGraphStorage ConstructStorage() const;
    };

}