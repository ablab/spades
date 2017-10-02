#pragma once
#include "common/assembly_graph/core/graph.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "scaffold_graph_extractor.hpp"
#include "common/barcode_index/cluster_storage.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"

namespace path_extend {
    class ScaffoldSubgraphExtractor {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef cluster_storage::Cluster::SimpleGraph<ScaffoldVertex> SimpleGraph;

        virtual SimpleGraph ExtractSubgraphBetweenVertices(const ScaffoldGraph& scaffold_graph, const ScaffoldVertex& first,
                                                           const ScaffoldVertex& second) const = 0;
    };

    class CloudScaffoldSubgraphExtractor: public ScaffoldSubgraphExtractor {
     public:
        using ScaffoldSubgraphExtractor::ScaffoldGraph;
        using ScaffoldSubgraphExtractor::ScaffoldVertex;
        using ScaffoldSubgraphExtractor::ScaffoldEdge;
        using ScaffoldSubgraphExtractor::SimpleGraph;
     private:
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_;
        const size_t distance_threshold_;
        const double share_threshold_;
        const size_t count_threshold_;
        const size_t small_length_threshold_;
        const size_t large_length_threshold_;
     public:

     public:
        CloudScaffoldSubgraphExtractor(const Graph& g_,
                                        const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                                        const size_t distance_threshold_,
                                        double share_threshold,
                                        const size_t count_threshold_,
                                        const size_t small_length_threshold,
                                        const size_t large_length_threshold);

        SimpleGraph ExtractSubgraphBetweenVertices(const ScaffoldGraph& scaffold_graph, const ScaffoldVertex& first,
                                                   const ScaffoldVertex& second) const override;

        bool CheckSubGraphVertex (const ScaffoldVertex& vertex, const ScaffoldVertex& first, const ScaffoldVertex& second) const;

        bool CheckSubgraphEdge (const ScaffoldEdge& edge, const ScaffoldVertex& first,
                                const ScaffoldVertex& second, const unordered_set<ScaffoldVertex>& subgraph_vertices) const;

        SimpleGraph RemoveDisconnectedVertices(const SimpleGraph& graph, const EdgeId& source, const EdgeId& sink) const ;

        DECL_LOGGER("ScaffoldGraphSubgraphExtractor");
    };

    class ReachabilityChecker {
     public:
        typedef debruijn_graph::EdgeId VertexT;
        typedef cluster_storage::Cluster::SimpleGraph<VertexT> SimpleGraph;
     private:
        std::unordered_set<VertexT> visited_;
        std::unordered_set<VertexT> passed_;
     protected:
        const SimpleGraph& graph_;

     public:
        ReachabilityChecker(const SimpleGraph& graph_);
        virtual ~ReachabilityChecker();
        void Run(const VertexT& start, const VertexT& target);
        unordered_set<VertexT> GetPassedVertices();

        DECL_LOGGER("ReachabilityChecker");
     private:
        virtual SimpleGraph::const_iterator GetBeginIterator(
            const VertexT& vertex) const = 0;
        virtual SimpleGraph::const_iterator GetEndIterator(const VertexT& vertex) const = 0;

        bool ProcessVertex(const VertexT& vertex, const VertexT& target);
    };

    class ForwardReachabilityChecker: public ReachabilityChecker {
        using ReachabilityChecker::graph_;
        using ReachabilityChecker::SimpleGraph;
     public:
        ForwardReachabilityChecker(const SimpleGraph& graph_);
     private:
        SimpleGraph::const_iterator GetBeginIterator(
            const VertexT& vertex) const override;
        SimpleGraph::const_iterator GetEndIterator(const VertexT& vertex) const override;
    };

    class BackwardReachabilityChecker: public ReachabilityChecker {
        using ReachabilityChecker::graph_;
        using ReachabilityChecker::SimpleGraph;
     public:
        BackwardReachabilityChecker(const SimpleGraph& graph_);
     private:
        SimpleGraph::const_iterator GetBeginIterator(const VertexT& vertex) const override;
        SimpleGraph::const_iterator GetEndIterator(const VertexT& vertex) const override;
    };

    class SubgraphPathExtractor {
     public:
        typedef debruijn_graph::EdgeId EdgeId;
        typedef cluster_storage::Cluster::SimpleGraph<EdgeId> SimpleGraph;

        vector<EdgeId> ExtractPathFromSubgraph(const SimpleGraph& graph, const EdgeId& source, const EdgeId& sink);

     private:

        vector<EdgeId> GetSimplePath(const SimpleGraph& graph, const EdgeId& source, const EdgeId& sink);
     public:
        DECL_LOGGER("SubgraphPathExtractor");
    };

    class InsertedVerticesData {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldVertex ScaffoldVertex;

        const std::unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_connections_map_;
        const size_t inserted_vertices_;
        const set<ScaffoldGraph::ScaffoldEdge> closed_edges_;
     public:
        const unordered_map<ScaffoldVertex, ScaffoldVertex>& GetInsertedConnectionsMap() const;
     public:
        size_t GetInsertedVertices() const;

        InsertedVerticesData(const unordered_map<ScaffoldVertex, ScaffoldVertex>& inserted_connections_map_,
                             const size_t inserted_vertices_, const std::set<ScaffoldGraph::ScaffoldEdge>& closed_edges);

        set<ScaffoldGraph::ScaffoldEdge> GetClosedEdges() const;
    };

    class IterationResult {
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef ScaffoldGraph::ScaffoldVertex ScaffoldVertex;

        const ScaffoldGraph new_graph_;
        const size_t inserted_vertices_;
        const std::set<ScaffoldEdge> closed_edges_;

     public:
        IterationResult(const ScaffoldGraph& new_graph_,
                        const size_t inserted_vertices_,
                        const std::set<ScaffoldEdge>& closed_edges_);

        const ScaffoldGraph& GetNewGraph() const;
        size_t GetInsertedVertices() const;
        std::set<ScaffoldEdge> GetClosedEdges() const;
    };

    class ScaffoldGraphGapCloser {
     public:
        typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
        typedef ScaffoldGraph::ScaffoldVertex ScaffoldVertex;
        typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
        typedef debruijn_graph::Graph Graph;

     private:
        const Graph& g_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t distance_threshold_;
        const double share_threshold_;
        const size_t count_threshold_;
        const size_t small_length_threshold_;
        const size_t large_length_threshold_;

     public:
        ScaffoldGraphGapCloser(const Graph& g_,
                               const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                               const size_t distance_threshold_,
                               const double share_threshold_,
                               const size_t count_threshold_,
                               const size_t small_length_threshold_,
                               const size_t large_length_threshold_);

        ScaffoldGraph ExtractGapClosingPaths(const ScaffoldGraph& large_scaffold_graph,
                                             const ScaffoldGraph& small_scaffold_graph) const;

        IterationResult LaunchGapClosingIteration(const ScaffoldGraph& current_graph,
                                                                   const vector<ScaffoldEdge>& univocal_edges) const;

        InsertedVerticesData GetInsertedConnections(const vector<ScaffoldEdge>& univocal_edges,
                                                    const ScaffoldGraph& current_graph) const;
        DECL_LOGGER("ScaffoldGraphGapCloser");
    };
}