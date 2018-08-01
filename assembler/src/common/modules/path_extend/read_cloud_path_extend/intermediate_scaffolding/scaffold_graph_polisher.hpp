#pragma once
#include "common/assembly_graph/core/graph.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/barcode_index/cluster_storage.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "pe_extraction.hpp"
#include "predicate_builders.hpp"
#include "path_cluster_helper.hpp"

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

    struct PathExtractionParams {
      const size_t linkage_distance_;
      const double path_cluster_relative_threshold_;
      const size_t min_read_threshold_;

      PathExtractionParams(size_t linkage_distance, double path_cluster_relative_threshold, size_t min_read_threshold);
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
        explicit ReachabilityChecker(const SimpleTransitionGraph& graph_);
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
        explicit ForwardReachabilityChecker(const SimpleTransitionGraph& graph_);
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

    class ScaffoldGraphGapCloserParamsConstructor {
     public:
        CloudSubgraphExtractorParams ConstructSubgraphExtractorParamsFromConfig(size_t length_upper_bound);
        PathExtractionParams ConstructPathExtractorParamsFromConfig();
    };

    class ScaffoldIndexInfoExtractorHelper {
     public:
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> ConstructIndexExtractorFromParams(
            const scaffold_graph::ScaffoldGraph scaffold_graph,
            const conj_graph_pack &gp,
            const CloudSubgraphExtractorParams& subgraph_extractor_params) const;
    };

    class ScaffoldGraphPolisher {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef debruijn_graph::Graph Graph;
        typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
        typedef vector<vector<ScaffoldVertex>> InternalPaths;

     private:
        const Graph &g_;
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor_;
        shared_ptr<CorrectPathExtractor> path_extractor_;
        const CloudSubgraphExtractorParams &subgraph_extractor_params_;

     public:
        ScaffoldGraphPolisher(const Graph &g_,
                              shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor,
                              shared_ptr<CorrectPathExtractor> path_extractor,
                              const CloudSubgraphExtractorParams &subgraph_extractor_params);

        ScaffoldGraph CleanSmallGraphUsingLargeGraph(const ScaffoldGraph &large_scaffold_graph,
                                                     const ScaffoldGraph &small_scaffold_graph) const;

        InternalPaths ExtractPathsWithinUnivocal(const ScaffoldGraph &current_graph,
                                                 const vector<ScaffoldEdge> &univocal_edges) const;

//        ScaffoldGraph CleanGraphUsingCutVertices(const ScaffoldGraph& input_graph, const vector<ScaffoldEdge>& univocal_edges) const;
        DECL_LOGGER("ScaffoldGraphPolisher");
    };

    class ScaffoldGraphPolisherLauncher {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

     public:
        ScaffoldGraph GetFinalScaffoldGraph(const conj_graph_pack &graph_pack,
                                            const ScaffoldGraphStorage &scaffold_graph_storage, bool path_scaffolding);

     private:
        shared_ptr<cluster_storage::InitialClusterStorage> ConstructInitialStorage(const conj_graph_pack &gp,
                                                                                   const ScaffoldGraph &scaffold_graph,
                                                                                   const PathExtractionParams &params,
                                                                                   bool path_scaffolding) const;
    };
}