#pragma once
#include "common/modules/path_extend/scaffolder2015/connection_condition2015.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"

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
        virtual bool IsLast() const override;
    };

    struct ScaffolderParams {
      size_t tail_threshold_;
      size_t count_threshold_;
      double score_threshold_;
      size_t barcode_threshold_;
      size_t length_threshold_;
      size_t initial_distance_;
      double middle_fraction_;
      double split_procedure_strictness_;
    };

    class ScaffoldEdgePredicate {
     public:
        typedef scaffold_graph::ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

        virtual bool Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const = 0;
        virtual ~ScaffoldEdgePredicate() = default;
    };

    struct LongGapDijkstraParams {
      const size_t barcode_threshold_;
      const size_t count_threshold_;
      const size_t tail_threshold_;
      const size_t len_threshold_;
      const size_t distance_;

      LongGapDijkstraParams(const size_t barcode_threshold_,
                                  const size_t count_threshold_,
                                  const size_t tail_threshold_,
                                  const size_t len_threshold_,
                                  const size_t distance);
    };

    class LongGapDijkstraPredicate: public ScaffoldEdgePredicate {
        using ScaffoldEdgePredicate::ScaffoldEdge;

        const Graph& g;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const LongGapDijkstraParams params_;
     public:
        LongGapDijkstraPredicate(const Graph& g,
                                         const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                         const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                         const LongGapDijkstraParams& params);
        bool Check(const ScaffoldEdge& scaffold_edge) const override;

    };

    class EdgeSplitPredicate: public ScaffoldEdgePredicate {
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef barcode_index::BarcodeId BarcodeId;

        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t count_threshold_;
        const double strictness_;
     public:
        EdgeSplitPredicate(const Graph& g_,
                           const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                           const size_t count_threshold_,
                           double strictness);

        bool Check(const ScaffoldEdge& scaffold_edge) const override;

     private:
        bool CheckOrderingForThreeSegments(const vector<BarcodeId>& first, const vector<BarcodeId>& second,
                                           const vector<BarcodeId>& third, double strictness) const;

        bool CheckOrderingForFourSegments(const vector<BarcodeId>& first, const vector<BarcodeId>& second,
                                          const vector<BarcodeId>& third, const vector<BarcodeId>& fourth) const;

        DECL_LOGGER("EdgeSplitPredicate");
    };

    class EdgeInTheMiddlePredicate {
     public:
        typedef barcode_index::BarcodeId BarcodeId;

     private:
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t count_threshold_;
        const double shared_fraction_threshold_;

     public:
        EdgeInTheMiddlePredicate(const Graph& g_,
                                         const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                         size_t count_threshold,
                                         double shared_fraction_threshold);

        bool IsCorrectOrdering(const EdgeId& first, const EdgeId& second, const EdgeId& third);
        DECL_LOGGER("EdgeInTheMiddlePredicate");
    };

    class NextFarEdgesPredicate: public ScaffoldEdgePredicate {
     public:
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef scaffold_graph::ScaffoldGraph::ScaffoldVertex ScaffoldVertex;

     private:
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t count_threshold_;
        const double shared_fraction_threshold_;
        const std::function<vector<ScaffoldVertex>(ScaffoldVertex)>& candidates_getter_;
     public:
        NextFarEdgesPredicate(const Graph& g_,
                              const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                              const size_t count_threshold_,
                              const double shared_fraction_threshold_,
                              const std::function<vector<ScaffoldVertex>(ScaffoldVertex)>& candidates_getter_);

        bool Check(const ScaffoldEdge& scaffold_edge) const override;
    };

    class TransitiveEdgesPredicate: public ScaffoldEdgePredicate {
     public:
        using ScaffoldEdgePredicate::ScaffoldEdge;
        typedef scaffold_graph::ScaffoldGraph::ScaffoldVertex ScaffoldVertex;

     private:
        const std::function<vector<ScaffoldVertex>(ScaffoldVertex)>& candidates_getter_;
     public:
        TransitiveEdgesPredicate(const std::function<vector<ScaffoldVertex>(ScaffoldVertex)>& candidates_getter_);

        bool Check(const ScaffoldEdge& scaffold_edge) const override;
    };

    class EdgePairScoreFunction {
     public:
        virtual double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const = 0;
        virtual ~EdgePairScoreFunction() = default;
    };

    class BarcodeScoreFunction: public EdgePairScoreFunction {
        const size_t read_count_threshold_;
        const size_t tail_threshold_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const Graph& graph_;

     public:
        BarcodeScoreFunction(const size_t read_count_threshold,
                             const size_t tail_threshold,
                             const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                             const Graph& graph);

        virtual double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const;
    };
}
