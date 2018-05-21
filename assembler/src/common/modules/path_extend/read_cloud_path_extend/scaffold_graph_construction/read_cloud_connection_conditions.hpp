#pragma once
#include "common/modules/path_extend/path_extender.hpp"
#include "common/modules/path_extend/pipeline/launch_support.hpp"
#include "common/modules/path_extend/extension_chooser.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

namespace path_extend {
    //Same as AssemblyGraphConnectionCondition, but stops after reaching unique edges.
    class AssemblyGraphUniqueConnectionCondition : public AssemblyGraphConnectionCondition {
        using AssemblyGraphConnectionCondition::g_;
        using AssemblyGraphConnectionCondition::interesting_edge_set_;
        using AssemblyGraphConnectionCondition::max_connection_length_;
        //fixme duplication with interesting edges, needed to pass to dijkstra
        const ScaffoldingUniqueEdgeStorage& unique_storage_;
     public:
        AssemblyGraphUniqueConnectionCondition(const Graph& g,
                                               size_t max_connection_length,
                                               const ScaffoldingUniqueEdgeStorage& unique_edges);
        map<EdgeId, double> ConnectedWith(EdgeId e) const override;
        bool IsLast() const override;
    };

    class ScaffoldEdgePredicate: public func::AbstractPredicate<const scaffold_graph::ScaffoldGraph::ScaffoldEdge&> {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

        virtual ~ScaffoldEdgePredicate() = default;
    };

    struct ReadCloudMiddleDijkstraParams {
      const size_t count_threshold_;
      const size_t tail_threshold_;
      const size_t distance_;

      const LongEdgePairGapCloserParams edge_pair_gap_closer_params_;

      ReadCloudMiddleDijkstraParams(size_t count_threshold_,
                                    size_t tail_threshold_,
                                    size_t distance_,
                                    const LongEdgePairGapCloserParams &edge_pair_gap_closer_params_);
    };

    class ReadCloudMiddleDijkstraPredicate: public ScaffoldEdgePredicate {
        using ScaffoldEdgePredicate::ScaffoldEdge;

        const Graph& g;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        const shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor_;
        const shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
        const ReadCloudMiddleDijkstraParams params_;
     public:
        ReadCloudMiddleDijkstraPredicate(const Graph& g,
                                         const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                         shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
                                         shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                                         const ReadCloudMiddleDijkstraParams& params);
        bool Check(const ScaffoldEdge& scaffold_edge) const override;

        DECL_LOGGER("ReadCloudMiddleDijkstraPredicate");
    };

    struct CompositeConnectionParams {
      const size_t paired_lib_index_;
      const size_t prefix_length_;
      const config::dataset& dataset_info;
      const path_extend::PathExtendParamsContainer pe_params_;

      CompositeConnectionParams(size_t paired_lib_index_,
                                size_t prefix_length_,
                                const config::dataset &dataset_info,
                                const PathExtendParamsContainer &pe_params_);
    };

    class CompositeConnectionPredicate: public ScaffoldEdgePredicate {
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;

        const conj_graph_pack& gp_;
        shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        const omnigraph::de::PairedInfoIndicesT<Graph>& clustered_indices_;
        const size_t length_bound_;
        const CompositeConnectionParams params_;
        const LongEdgePairGapCloserParams predicate_params_;
        bool scaffolding_mode_;

     public:

        CompositeConnectionPredicate(const conj_graph_pack &gp_,
                                     shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> short_edge_extractor,
                                     shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                     const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                     const de::PairedInfoIndicesT<debruijn_graph::DeBruijnGraph> &clustered_indices_,
                                     size_t length_bound_,
                                     const CompositeConnectionParams &params_,
                                     const LongEdgePairGapCloserParams &predicate_params_,
                                     bool scaffolding_mode_);

        bool Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const override;

     private:
        bool SearchTargetUsingExtenders(QueueContainer& paths_container,
                                        GraphCoverageMap& cover_map,
                                        vector<shared_ptr<SearchingMultiExtender>> extenders,
                                        const VertexId& target_vertex) const;

        //fixme huge code duplication with ExtendersGenerator
        shared_ptr<path_extend::ExtensionChooser> ConstructSimpleExtensionChooser() const;
        shared_ptr<path_extend::ExtensionChooser> ConstructScaffoldingExtensionChooser() const;

        shared_ptr<ScaffoldVertexPredicate> ConstructScaffoldVertexPredicate(
            const ScaffoldVertex& start, const ScaffoldVertex& end,
            shared_ptr<path_extend::PairEntryProcessor> entry_processor) const;

        shared_ptr<path_extend::SearchingMultiExtender> ConstructBasicSearchingExtender(
            path_extend::QueueContainer& paths_container,
            GraphCoverageMap& cover_map,
            shared_ptr<ExtensionChooser> extension_chooser) const;

        shared_ptr<path_extend::SearchingMultiExtender> ConstructScaffoldingSearchingExtender(
            path_extend::QueueContainer& paths_container,
            GraphCoverageMap& cover_map,
            shared_ptr<ExtensionChooser> extension_chooser) const;

        shared_ptr<GapAnalyzer> MakeGapAnalyzer(double is_variation) const;

        DECL_LOGGER("CompositeConnectionPredicate");
    };

    class EdgeSplitPredicate: public ScaffoldEdgePredicate {
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef barcode_index::BarcodeId BarcodeId;

        const Graph& g_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_;
        const size_t count_threshold_;
        const double strictness_;
     public:
        EdgeSplitPredicate(const Graph& g_,
                           shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
                           size_t count_threshold_,
                           double strictness);

        bool Check(const ScaffoldEdge& scaffold_edge) const override;

     private:
        bool CheckOrderingForThreeSegments(const barcode_index::SimpleVertexEntry& first,
                                           const barcode_index::SimpleVertexEntry& second,
                                           const barcode_index::SimpleVertexEntry& third, double strictness) const;

        bool CheckOrderingForFourSegments(const barcode_index::SimpleVertexEntry &first,
                                          const barcode_index::SimpleVertexEntry &second,
                                          const barcode_index::SimpleVertexEntry &third,
                                          const barcode_index::SimpleVertexEntry &fourth) const;

        DECL_LOGGER("EdgeSplitPredicate");
    };

    class EdgeInTheMiddlePredicate {
     public:
        typedef barcode_index::BarcodeId BarcodeId;

     private:
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_;
        const size_t count_threshold_;
        const double shared_fraction_threshold_;

     public:
        EdgeInTheMiddlePredicate(const Graph& g_,
                                 const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                 size_t count_threshold,
                                 double shared_fraction_threshold);

        bool IsCorrectOrdering(const EdgeId& first, const EdgeId& second, const EdgeId& third);
        DECL_LOGGER("EdgeInTheMiddlePredicate");
    };

    class SimpleSearcher {
     public:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
     private:
        const ScaffoldGraph& scaff_graph_;
        const Graph& g_;
        size_t distance_threshold_;

        struct VertexWithDistance {
          ScaffoldVertex vertex;
          size_t distance;
          VertexWithDistance(const ScaffoldVertex& vertex, size_t distance);
        };

     public:
        SimpleSearcher(const scaffold_graph::ScaffoldGraph& graph_, const Graph& g, size_t distance_);

        vector<ScaffoldVertex> GetReachableVertices(const ScaffoldVertex& vertex, const ScaffoldGraph::ScaffoldEdge& restricted_edge);

        void ProcessVertex(std::queue<VertexWithDistance>& vertex_queue, const VertexWithDistance& vertex,
                           std::unordered_set<ScaffoldVertex>& visited, const ScaffoldGraph::ScaffoldEdge& restricted_edge);

        bool AreEqual(const ScaffoldGraph::ScaffoldEdge& first, const ScaffoldGraph::ScaffoldEdge& second);

        DECL_LOGGER("SimpleSearcher");
    };

    class TransitiveEdgesPredicate: public ScaffoldEdgePredicate {
     public:
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     private:
        const scaffold_graph::ScaffoldGraph scaffold_graph_;
        const Graph& g_;
        size_t distance_threshold_;
     public:
        TransitiveEdgesPredicate(const scaffold_graph::ScaffoldGraph& graph, const Graph& g, size_t distance_threshold);

        bool Check(const ScaffoldEdge& scaffold_edge) const override;

        DECL_LOGGER("TransitiveEdgesPredicate");
    };

    class ScaffoldEdgeScoreFunction {
     protected:
        typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
     public:
        virtual double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const = 0;
        virtual ~ScaffoldEdgeScoreFunction() = default;
    };

    class AbstractBarcodeScoreFunction: public ScaffoldEdgeScoreFunction {
     protected:
        const Graph& graph_;
        shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_;

     public:
        AbstractBarcodeScoreFunction(
            const Graph& graph_,
            shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_);
    };

    class NormalizedBarcodeScoreFunction: public AbstractBarcodeScoreFunction {
        using AbstractBarcodeScoreFunction::barcode_extractor_;
        using AbstractBarcodeScoreFunction::graph_;

     public:
        NormalizedBarcodeScoreFunction(const Graph& graph_,
                                       shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_);

        double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const override;

        DECL_LOGGER("NormalizedBarcodeScoreFunction");
    };

    class TrivialBarcodeScoreFunction: public AbstractBarcodeScoreFunction {
        using AbstractBarcodeScoreFunction::barcode_extractor_;
        using AbstractBarcodeScoreFunction::graph_;
        const size_t read_count_threshold_;
        const size_t tail_threshold_;

     public:
        TrivialBarcodeScoreFunction(
            const Graph& graph_,
            shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
            size_t read_count_threshold_,
            size_t tail_threshold_);

        double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const override;
    };

}
