#pragma once
#include "common/assembly_graph/core/graph.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/barcode_index/cluster_storage.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "pe_extraction.hpp"
#include "predicate_builders.hpp"

namespace path_extend {
    struct CloudSubgraphExtractorParams {
      const size_t distance_threshold_;
      const double share_threshold_;
      const size_t count_threshold_;
      const size_t small_length_threshold_;
      const size_t large_length_threshold_;

      CloudSubgraphExtractorParams(
          size_t distance_threshold_,
          double share_threshold_,
          size_t count_threshold_,
          size_t small_length_threshold_,
          size_t large_length_threshold_);
    };

    struct PathClusterPredicateParams {
      const size_t linkage_distance_;
      const double path_cluster_threshold_;
      const size_t min_read_threshold_;

      PathClusterPredicateParams(size_t linkage_distance_,
                                 double path_cluster_threshold_,
                                 size_t min_read_threshold_);
    };

    struct PathExtractorParts {
      const vector<shared_ptr<GapCloserPredicateBuilder>>& predicate_builders_;
      const shared_ptr<GapCloserScoreFunctionBuilder> score_builder_;

      PathExtractorParts(const vector<shared_ptr<GapCloserPredicateBuilder>>& predicate_builders_,
                          const shared_ptr<GapCloserScoreFunctionBuilder>& score_builder_);
    };

    class ScaffoldSubgraphExtractor {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::VertexId ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef path_extend::SimpleGraph<ScaffoldVertex> SimpleGraph;

        virtual SimpleGraph ExtractSubgraphBetweenVertices(const ScaffoldGraph& scaffold_graph, const ScaffoldVertex& first,
                                                           const ScaffoldVertex& second) const = 0;
    };

    class GapCloserUtils {
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef path_extend::GapCloserPredicateBuilder PredicateBuilder;
        typedef SimpleGraph<path_extend::scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;
     public:
        bool IsSimplePath(const SimpleTransitionGraph& graph, const ScaffoldVertex& source, const ScaffoldVertex& sink) const;

        SimpleTransitionGraph RemoveDisconnectedVertices(const SimpleTransitionGraph& graph, const ScaffoldVertex& source,
                                                         const ScaffoldVertex& sink) const;

//        std::unordered_set<ScaffoldVertex> ExtractCutVertices(const SimpleTransitionGraph& graph,
//                                                              const ScaffoldVertex& source,
//                                                              const ScaffoldVertex& sink) const;
    };

    class CutVerticesExtractor {
     public:
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

     private:
        const SimpleTransitionGraph& graph_;
     public:
        CutVerticesExtractor(const SimpleTransitionGraph &graph_);

        vector<ScaffoldVertex> GetCutVertices(const ScaffoldVertex& source, const ScaffoldVertex& sink);

     private:
        bool Check(const ScaffoldVertex& sink, const ScaffoldVertex& source, const ScaffoldVertex& candidate);

        DECL_LOGGER("CutVerticesExtractor");
    };

    class CloudScaffoldSubgraphExtractor: public ScaffoldSubgraphExtractor {
     public:
        using ScaffoldSubgraphExtractor::ScaffoldGraph;
        using ScaffoldSubgraphExtractor::ScaffoldVertex;
        using ScaffoldSubgraphExtractor::ScaffoldEdge;
        using ScaffoldSubgraphExtractor::SimpleGraph;
     private:
        const Graph& g_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor_;
        const CloudSubgraphExtractorParams params_;
     public:

     public:
        CloudScaffoldSubgraphExtractor(const Graph& g_,
                                       shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> extractor_,
                                       const CloudSubgraphExtractorParams& params);

        SimpleGraph ExtractSubgraphBetweenVertices(const ScaffoldGraph& scaffold_graph, const ScaffoldVertex& first,
                                                   const ScaffoldVertex& second) const override;

        bool CheckSubGraphVertex (const ScaffoldVertex& vertex, const ScaffoldVertex& first, const ScaffoldVertex& second) const;

        bool CheckSubgraphEdge (const ScaffoldEdge& edge, const ScaffoldVertex& first,
                                const ScaffoldVertex& second, const unordered_set<ScaffoldVertex>& subgraph_vertices) const;

        DECL_LOGGER("ScaffoldGraphSubgraphExtractor");
    };

    class ReachabilityChecker {
     public:
        typedef path_extend::scaffold_graph::ScaffoldVertex VertexT;
        typedef path_extend::SimpleGraph<VertexT> SimpleTransitionGraph;
     private:
        std::unordered_set<VertexT> visited_;
        std::unordered_set<VertexT> passed_;
     protected:
        const SimpleTransitionGraph& graph_;

     public:
        ReachabilityChecker(const SimpleTransitionGraph& graph_);
        virtual ~ReachabilityChecker();
        void Run(const VertexT& start, const VertexT& target);
        unordered_set<VertexT> GetPassedVertices();

        DECL_LOGGER("ReachabilityChecker");
     private:
        virtual SimpleTransitionGraph::const_iterator GetBeginIterator(
            const VertexT& vertex) const = 0;
        virtual SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT& vertex) const = 0;

        bool ProcessVertex(const VertexT& vertex, const VertexT& target);
    };

    class ForwardReachabilityChecker: public ReachabilityChecker {
        using ReachabilityChecker::graph_;
        using ReachabilityChecker::SimpleTransitionGraph;
     public:
        ForwardReachabilityChecker(const SimpleTransitionGraph& graph_);
     private:
        SimpleTransitionGraph::const_iterator GetBeginIterator(
            const VertexT& vertex) const override;
        SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT& vertex) const override;
    };

    class BackwardReachabilityChecker: public ReachabilityChecker {
        using ReachabilityChecker::graph_;
        using ReachabilityChecker::SimpleTransitionGraph;
     public:
        BackwardReachabilityChecker(const SimpleTransitionGraph& graph_);
     private:
        SimpleTransitionGraph::const_iterator GetBeginIterator(const VertexT& vertex) const override;
        SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT& vertex) const override;
    };

    class SubgraphPathExtractor {
     public:
        typedef vector<shared_ptr<GapCloserPredicateBuilder>> p_builders_t;
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
     private:
        vector<shared_ptr<GapCloserPredicateBuilder>> predicate_builders_;
        //todo multiple scores?
        shared_ptr<GapCloserScoreFunctionBuilder> score_function_builder_;

     public:
        typedef SimpleGraph<path_extend::scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

        SubgraphPathExtractor(const p_builders_t& predicate_builders,
                              shared_ptr<GapCloserScoreFunctionBuilder> score_function_builder);

        vector<ScaffoldVertex> ExtractPathFromSubgraph(const SimpleTransitionGraph& graph,
                                                       const ScaffoldVertex& source,
                                                       const ScaffoldVertex& sink) const;

        vector<ScaffoldVertex> ExtractSimplePathFromSubgraph(const SimpleTransitionGraph& graph,
                                                             const ScaffoldVertex& source,
                                                             const ScaffoldVertex& sink) const;

        vector<ScaffoldVertex> ExtractPathUsingScoreFunction(const SimpleTransitionGraph& graph, const ScaffoldVertex& source,
                                                             const ScaffoldVertex& sink,
                                                             shared_ptr<ScaffoldEdgeScoreFunction> score_function) const;

     private:

        std::pair<ScaffoldVertex, double> GetNextMaxEdge(const ScaffoldVertex& current,
                                                         shared_ptr<ScaffoldEdgeScoreFunction> score_function,
                                                         const SubgraphPathExtractor::SimpleTransitionGraph& graph) const;

        std::pair<ScaffoldVertex, double> GetPrevMaxEdge(const ScaffoldVertex& current,
                                                         shared_ptr<ScaffoldEdgeScoreFunction> score_function,
                                                         const SubgraphPathExtractor::SimpleTransitionGraph& graph) const;

        vector<ScaffoldVertex> GetSimplePath(const SimpleTransitionGraph& graph, const ScaffoldVertex& source, const ScaffoldVertex& sink) const;

     public:
        DECL_LOGGER("SubgraphPathExtractor");
    };

    class SubgraphEdgeChecker {
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
        typedef vector<shared_ptr<GapCloserPredicateBuilder>> p_builders_t;


     public:
        SimpleTransitionGraph CleanGraphUsingPredicate(SimpleTransitionGraph& graph,
                                                       shared_ptr<ScaffoldEdgePredicate> predicate_ptr) const;

        SimpleTransitionGraph CleanGraphUsingPredicateBuilders(SimpleTransitionGraph& graph, const ScaffoldVertex& source,
                                                               const ScaffoldVertex& sink,
                                                               const p_builders_t& predicate_builders) const;
    };

    class InsertedVerticesData {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

        const std::unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_connections_map_;
        const size_t inserted_vertices_;
        const set<ScaffoldGraph::ScaffoldEdge> closed_edges_;
     public:
        const unordered_map<ScaffoldVertex, ScaffoldVertex>& GetInsertedConnectionsMap() const;
     public:
        size_t GetInsertedVertices() const;

        InsertedVerticesData(const unordered_map<ScaffoldVertex, ScaffoldVertex>& inserted_connections_map_,
                             size_t inserted_vertices_, const std::set<ScaffoldGraph::ScaffoldEdge>& closed_edges);

        set<ScaffoldGraph::ScaffoldEdge> GetClosedEdges() const;
    };

    class IterationResult {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

        const ScaffoldGraph new_graph_;
        const size_t inserted_vertices_;
        const std::set<ScaffoldEdge> closed_edges_;

     public:
        IterationResult(const ScaffoldGraph& new_graph_,
                        size_t inserted_vertices_,
                        const std::set<ScaffoldEdge>& closed_edges_);

        const ScaffoldGraph& GetNewGraph() const;
        size_t GetInsertedVertices() const;
        std::set<ScaffoldEdge> GetClosedEdges() const;
    };

    class ScaffoldGraphGapCloserParamsConstructor {
     public:
        CloudSubgraphExtractorParams ConstructSubgraphExtractorParamsFromConfig();
        PathClusterPredicateParams ConstructPathClusterPredicateParamsFromConfig();
    };

    class PathExtractionPartsConstructor {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        const conj_graph_pack& gp_;
     public:
        explicit PathExtractionPartsConstructor(const conj_graph_pack& gp_);
        vector<shared_ptr<GapCloserPredicateBuilder>> ConstructPredicateBuilders() const;
        shared_ptr<GapCloserScoreFunctionBuilder> ConstructPathClusterScoreFunction(const PathClusterPredicateParams& params) const;
        shared_ptr<GapCloserPredicateBuilder> ConstructPEPredicate() const;
    };

    class ScaffoldGraphGapCloser {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef debruijn_graph::Graph Graph;
        typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

     private:
        const Graph& g_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor_;
        const CloudSubgraphExtractorParams& subgraph_extractor_params_;
        const PathExtractorParts& path_extractor_params_;

     public:
        ScaffoldGraphGapCloser(const Graph& g_,
                               shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor,
                               const CloudSubgraphExtractorParams& subgraph_extractor_params,
                               const PathExtractorParts& path_extractor_params);

        ScaffoldGraph CleanSmallGraphUsingLargeGraph(const ScaffoldGraph &large_scaffold_graph,
                                                     const ScaffoldGraph &small_scaffold_graph) const;

        ScaffoldGraph CloseGapsInLargeGraph(const ScaffoldGraph& large_scaffold_graph,
                                            const ScaffoldGraph& small_scaffold_graph) const;

        IterationResult LaunchGapClosingIteration(const ScaffoldGraph& current_graph,
                                                                   const vector<ScaffoldEdge>& univocal_edges) const;

        InsertedVerticesData GetInsertedConnections(const vector<ScaffoldEdge>& univocal_edges,
                                                    const ScaffoldGraph& current_graph) const;

//        ScaffoldGraph CleanGraphUsingCutVertices(const ScaffoldGraph& input_graph, const vector<ScaffoldEdge>& univocal_edges) const;
        DECL_LOGGER("ScaffoldGraphGapCloser");
    };

    class ScaffoldGraphGapCloserLauncher {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

     public:
        ScaffoldGraph GetFinalScaffoldGraph(const conj_graph_pack& graph_pack);
    };
}