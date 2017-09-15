#include "common/barcode_index/cluster_storage_extractor.hpp"

namespace read_cloud_statistics {
class ClusterGraphAnalyzerTester {
 public:
    typedef cluster_storage::Cluster::InternalGraph InternalGraph;

    struct TestableGraph {
      InternalGraph graph_;
      string name_;

      TestableGraph(const InternalGraph& graph_, const string& name_) : graph_(graph_), name_(name_) {}
    };
 private:
    const debruijn_graph::conj_graph_pack& gp_;
    const contracted_graph::ContractedGraphFactoryHelper contracted_graph_builder_;

 public:
    ClusterGraphAnalyzerTester(const debruijn_graph::conj_graph_pack& gp_,
                               const contracted_graph::ContractedGraphFactoryHelper& contracted_graph_builder) :
        gp_(gp_), contracted_graph_builder_(contracted_graph_builder) {}

    //fixme move to unit tests
    void LaunchTests() {
        vector<TestableGraph> tests;
        tests.push_back(GenerateLinearGraph());
        tests.push_back(GenerateCycle());
        tests.push_back(GenerateSelfLoop());
        tests.push_back(GenerateSelfLoopWithEdge());
        tests.push_back(GenerateOneLoopPath());
        tests.push_back(GenerateTwoLoopPath());
        auto cluster_graph_analyzer = cluster_storage::ClusterGraphAnalyzer(contracted_graph_builder_);
        for (const auto& test: tests) {
            INFO("TESTING " << test.name_);
            TestGraph(test.graph_, cluster_graph_analyzer);
        }
    }

    vector<EdgeId> GetEdges(size_t number_of_edges) {
        unordered_set<EdgeId> edges;
        unordered_set<VertexId> vertices;
        omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(gp_.g);
        size_t counter = 0;
        for (const auto& edge: edge_iteration_helper) {
            if (CheckEdge(edge, vertices, edges)) {
                edges.insert(edge);
                vertices.insert(gp_.g.EdgeEnd(edge));
                vertices.insert(gp_.g.EdgeStart(edge));
                ++counter;
            }
            if (counter >= number_of_edges) {
                break;
            }
        }
        VERIFY(edges.size() * 2 == vertices.size());
        vector<EdgeId> result;
        std::copy(edges.begin(), edges.end(), std::back_inserter(result));
        return result;
    }

    bool CheckEdge(const EdgeId& edge, const unordered_set<VertexId>& vertices, const unordered_set<EdgeId>& edges) {
        bool conjugate_not_visited = edges.find(gp_.g.conjugate(edge)) == edges.end();
        bool not_loop = gp_.g.EdgeEnd(edge) != gp_.g.EdgeStart(edge);
        bool start_not_visited = vertices.find(gp_.g.EdgeStart(edge)) == vertices.end();
        bool end_not_visited = vertices.find(gp_.g.EdgeEnd(edge)) == vertices.end();
        return conjugate_not_visited and not_loop and start_not_visited and end_not_visited;
    }

    TestableGraph GenerateLinearGraph() {
        size_t number_of_edges = 4;
        auto edges = GetEdges(number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        for (auto first = edges.begin(), second = std::next(edges.begin()); second != edges.end(); ++first, ++second) {
            EdgeId start = *first;
            EdgeId end = *second;
            graph.AddEdge(start, end);
        }
        string name = "Linear graph";
        TestableGraph result(graph, name);
        return result;
    }

    TestableGraph GenerateCycle() {
        size_t number_of_edges = 4;
        auto edges = GetEdges(number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        for (auto first = edges.begin(), second = std::next(edges.begin()); second != edges.end(); ++first, ++second) {
            EdgeId start = *first;
            EdgeId end = *second;
            graph.AddEdge(start, end);
        }
        EdgeId start = edges.back();
        EdgeId end = edges[0];

        string name = "Cycle";
        graph.AddEdge(start, end);
        TestableGraph result(graph, name);
        return result;
    }

    TestableGraph GenerateSelfLoop() {
        size_t number_of_edges = 1;
        auto edges = GetEdges(number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        graph.AddEdge(edges[0], edges[0]);

        string name = "Self loop";
        TestableGraph result(graph, name);
        return result;
    }

    TestableGraph GenerateSelfLoopWithEdge() {
        size_t number_of_edges = 2;
        auto edges = GetEdges(number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        graph.AddEdge(edges[0], edges[0]);
        graph.AddEdge(edges[0], edges[1]);

        string name = "Self loop with edge";
        TestableGraph result(graph, name);
        return result;
    }

    TestableGraph GenerateOneLoopPath() {
        size_t number_of_edges = 4;
        auto edges = GetEdges(number_of_edges);
        VERIFY(edges.size() == number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        graph.AddEdge(edges[0], edges[1]);
        graph.AddEdge(edges[0], edges[3]);
        graph.AddEdge(edges[1], edges[2]);
        graph.AddEdge(edges[2], edges[3]);
        graph.AddEdge(edges[2], edges[1]);
        string name = "One loop path";
        TestableGraph result(graph, name);
        return result;
    }

    TestableGraph GenerateTwoLoopPath() {
        size_t number_of_edges = 6;
        auto edges = GetEdges(number_of_edges);
        VERIFY(edges.size() == number_of_edges);
        InternalGraph graph;
        for (const auto& edge: edges) {
            graph.AddVertex(edge);
        }
        graph.AddEdge(edges[0], edges[1]);
        graph.AddEdge(edges[0], edges[3]);
        graph.AddEdge(edges[0], edges[5]);
        graph.AddEdge(edges[1], edges[2]);
        graph.AddEdge(edges[2], edges[1]);
        graph.AddEdge(edges[2], edges[3]);
        graph.AddEdge(edges[2], edges[5]);
        graph.AddEdge(edges[3], edges[4]);
        graph.AddEdge(edges[4], edges[1]);
        graph.AddEdge(edges[4], edges[3]);
        graph.AddEdge(edges[4], edges[5]);
        string name = "Two loop path";
        TestableGraph result(graph, name);
        return result;
    }

    void TestGraph(const InternalGraph& graph, const cluster_storage::ClusterGraphAnalyzer analyzer) {
        auto contracted_graph = contracted_graph_builder_.ConstructFromInternalGraph(graph);
        cluster_storage::Cluster test_cluster(graph);
        string is_path = analyzer.IsPathCluster(test_cluster) ? "True" : "False";
        INFO("Is path cluster: " << is_path);
    }
};
}