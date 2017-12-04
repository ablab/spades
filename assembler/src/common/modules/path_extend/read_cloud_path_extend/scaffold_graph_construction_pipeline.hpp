#pragma once

#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_vertex.hpp"
#include "read_cloud_connection_conditions.hpp"
#include "scaffold_graph_storage.hpp"

namespace path_extend {
    struct ScaffolderParams {
      const size_t length_threshold_;
      size_t tail_threshold_;
      size_t count_threshold_;
      const double score_threshold_;
      const double connection_score_threshold_;
      const size_t connection_length_threshold_;
      const size_t connection_count_threshold_;
      size_t initial_distance_;
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
        ScaffolderParams ConstructScaffolderParamsFromCfg(size_t length_threshold) const;

        LongEdgePairGapCloserParams ConstructGapCloserParamsFromCfg(bool normalize_using_cov) const;
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
        shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;
        size_t max_threads_;

     public:
        BarcodeScoreConstructorCaller(const Graph& g_,
                                      shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                      size_t max_threads_);
        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class BarcodeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
        size_t max_threads;

     public:
        BarcodeConnectionConstructorCaller(const Graph& g_,
                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                           shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                           const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                           size_t max_threads);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class CompositeConnectionConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const conj_graph_pack& gp_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
        const size_t max_threads_;

     public:
        CompositeConnectionConstructorCaller(const conj_graph_pack &gp_,
                                             shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                             shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                             const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                             const size_t max_threads_);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class EdgeSplitConstructorCaller: public IterativeScaffoldGraphConstructorCaller {
        using IterativeScaffoldGraphConstructorCaller::ScaffoldGraph;
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
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
        size_t max_threads_;

     public:
        TransitiveConstructorCaller(const Graph& g_,
                                    size_t max_threads_);

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                        const ScaffoldGraph& scaffold_graph) const override;
    };

    class CloudScaffoldGraphConstructionPipeline {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> initial_constructor_;
        vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> construction_stages_;
        vector<shared_ptr<scaffold_graph::ScaffoldGraph>> intermediate_results_;
        const Graph &g_;
        const ScaffolderParams& params_;
        const string name_;

     public:
        CloudScaffoldGraphConstructionPipeline(shared_ptr<scaffold_graph::ScaffoldGraphConstructor> initial_constructor_,
                                               const Graph& g, const ScaffolderParams &params, const string &name);

        void AddStage(shared_ptr<IterativeScaffoldGraphConstructorCaller> stage);

        void Run();

        shared_ptr<scaffold_graph::ScaffoldGraph> GetResult() const;
    };

    class CloudScaffoldGraphConstuctor {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef scaffold_graph::ScaffoldGraphConstructor ScaffoldGraphConstructor;
        typedef path_extend::ScaffoldingUniqueEdgeStorage ScaffoldingUniqueEdgeStorage;
     private:
        const size_t max_threads_;
        const conj_graph_pack& gp_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

     public:
        CloudScaffoldGraphConstuctor(size_t max_threads_,
                                     const conj_graph_pack& gp,
                                     shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);
        ScaffoldGraph ConstructScaffoldGraphFromMinLength(size_t min_length) const;

        ScaffoldGraph ConstructScaffoldGraphFromMinLengthAndGraph(size_t min_length, const ScaffoldGraph& previous_graph) const;

        ScaffoldGraph ConstructScaffoldGraphFromPathContainer(const PathContainer& paths,
                                                              const ScaffoldingUniqueEdgeStorage& unique_storage,
                                                              size_t min_length) const;

        //todo replace storage with predicate
        ScaffoldGraph ConstructScaffoldGraphFromStorage(const ScaffolderParams& params,
                                                        const ScaffoldingUniqueEdgeStorage& unique_storage,
                                                        const set<ScaffoldVertex>& scaffold_vertices,
                                                        const string &initial_graph_name,
                                                        bool launch_full_pipeline,
                                                        bool path_merge_pipeline = false) const;

        ScaffoldGraph ConstructScaffoldGraphFromStorageAndGraph(const ScaffolderParams& params,
                                                                const ScaffoldGraph& previous_graph,
                                                                const ScaffoldingUniqueEdgeStorage& unique_storage,
                                                                const set<ScaffoldVertex>& scaffold_vertices,
                                                                bool launch_full_pipeline,
                                                                bool path_merge_pipeline = false) const;

     private:
        vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> ConstructStages(const ScaffolderParams& params,
                                                                                    const ScaffoldingUniqueEdgeStorage& unique_storage,
                                                                                    const set<ScaffoldVertex>& scaffold_vertices,
                                                                                    bool launch_full_pipeline,
                                                                                    bool path_merge_pipeline) const;
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