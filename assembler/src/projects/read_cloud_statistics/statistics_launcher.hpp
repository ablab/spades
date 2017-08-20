#pragma once

#include "common/barcode_index/contracted_graph.hpp"
#include "contracted_graph_stats/contracted_graph_analyzer.hpp"
#include "cluster_storage_analyzer.hpp"
#include "scaffold_graph.hpp"
#include "statistics_processor.hpp"
#include "transitions.hpp"
#include "../../common/modules/path_extend/pipeline/launch_support.hpp"
#include "../../common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "contig_path_analyzer.hpp"
#include "contracted_graph_stats/contracted_graph_local_statistics.hpp"

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
        const contracted_graph::ContractedGraphBuilder contracted_graph_builder_;

     public:
        ClusterGraphAnalyzerTester(const debruijn_graph::conj_graph_pack& gp_,
                                   const contracted_graph::ContractedGraphBuilder& contracted_graph_builder) :
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
            auto contracted_graph = contracted_graph_builder_.BuildContractedGraphFromInternalGraph(graph);
            cluster_storage::Cluster test_cluster(graph);
            string is_path = analyzer.IsPathCluster(test_cluster) ? "True" : "False";
            INFO("Is path cluster: " << is_path);
        }
    };

    struct StatisticsConfigs {
      const size_t edge_read_cluster_threshold;
      const size_t global_read_cluster_threshold;

      StatisticsConfigs(const size_t edge_read_cluster_threshold, const size_t global_read_cluster_threshold)
          : edge_read_cluster_threshold(edge_read_cluster_threshold),
            global_read_cluster_threshold(global_read_cluster_threshold) {}
    };

    class StatisticsLauncher {
        const debruijn_graph::conj_graph_pack& gp_;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
        StatisticsConfigs configs_;
     public:
        StatisticsLauncher(const debruijn_graph::conj_graph_pack& gp,
                           const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_,
                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                           const StatisticsConfigs &configs_) : gp_(gp), unique_storage_(unique_storage_),
                                                                barcode_extractor_ptr_(barcode_extractor_ptr), configs_(configs_) {}

        void Launch(const string& base_output_path) {
            INFO("Transition stats:");
            size_t distance = cfg::get().ts_res.distance;
            vector<size_t> distances = {1000, 2500, 5000, 7500, 10000, 20000, 35000, 50000};

            std::function<void(const string&, size_t)> path_cluster_analyzer = [=](const string& path, size_t dist) {
              this->AnalyzePathClusters(path, dist);
            };

            std::function<void(const string&, size_t)> contracted_graph_analyzer = [=](const string& path, size_t dist) {
              this->AnalyzeContractedGraph(path, dist);
            };


            INFO("Analyzing contig paths");
            AnalyzeContigPaths(base_output_path);

//            INFO("Analyzing path clusters");
//            LaunchAnalyzerForMultipleDistances(base_output_path, path_cluster_analyzer, distances, cfg::get().max_threads);
//            AnalyzePathClusters(base_output_path, distance);
//
//            AnalyzeTransitions(base_output_path, distance);
//            AnalyzeTransitionsForMultipleDistances(graph_pack, stats_path);
        }

     private:

        transitions::ContigTransitionStorage BuildFilteredTransitionStorage(const transitions::ContigPathBuilder& contig_path_builder,
                                                                            shared_ptr<transitions::TransitionStorageBuilder> transition_builder,
                                                                            const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                                                            const string& path_to_contigs) {
            const string EMPTY_STRING = "";
            transitions::ContigTransitionStorage result;
            if (path_to_contigs != EMPTY_STRING) {
                auto reference_paths = contig_path_builder.GetContigPaths(path_to_contigs);
                transitions::ContigPathFilter contig_path_filter(unique_storage);
                result = transition_builder->GetTransitionStorage(contig_path_filter.FilterPathsUsingUniqueStorage(
                    reference_paths));
            }
            INFO("Transition storage size: " << result.Size());
            return result;
        }

        std::unordered_map<string, transitions::ContigTransitionStorage>
        BuildTransitionStorages(const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage) {
            const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
            const string path_to_base_contigs = cfg::get().ts_res.statistics.base_contigs_path;
            const string path_to_cloud_contigs = cfg::get().ts_res.statistics.cloud_contigs_path;
            const string EMPTY_STRING = "";
            std::unordered_map<string, transitions::ContigTransitionStorage> name_to_transition_storage;
            INFO("Reference path: " << path_to_reference);
            INFO("Base contigs path: " << path_to_base_contigs);
            INFO("Cloud contigs path: " << path_to_cloud_contigs);

            transitions::ContigPathBuilder contig_path_builder(gp_);

            auto strict_transition_builder = make_shared<transitions::StrictTransitionStorageBuilder>(contig_path_builder);
            auto approximate_transition_builder = make_shared<transitions::ApproximateTransitionStorageBuilder>(contig_path_builder);

            auto reference_transition_storage =
                BuildFilteredTransitionStorage(contig_path_builder, strict_transition_builder,
                                               unique_storage, path_to_reference);
            name_to_transition_storage.insert({"Reference", reference_transition_storage});

            auto contig_transition_storage =
                BuildFilteredTransitionStorage(contig_path_builder, approximate_transition_builder,
                                               unique_storage, path_to_base_contigs);
            name_to_transition_storage.insert({"Contig", contig_transition_storage});

            auto cloud_transition_storage =
                BuildFilteredTransitionStorage(contig_path_builder, approximate_transition_builder,
                                               unique_storage, path_to_cloud_contigs);
            name_to_transition_storage.insert({"Read cloud contig", contig_transition_storage});

            return name_to_transition_storage;
        }

        void CountBasicClusterStats(const cluster_storage::ClusterStorage& cluster_storage,
                                    const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer) {
            size_t clusters = 0;
            size_t eulerian_clusters = 0;
            size_t path_clusters = 0;

            for (const auto& entry: cluster_storage) {
                auto cluster = entry.second;
                if (cluster.Size() >= 2 and cluster.GetReads() >= configs_.global_read_cluster_threshold) {
                    ++clusters;
                    if (ordering_analyzer.IsEulerianCluster(cluster)) {
                        ++eulerian_clusters;
                        if (ordering_analyzer.IsPathCluster(cluster)) {
                            TRACE("Is path cluster");
                            ++path_clusters;
                        } else {
                            TRACE("Not a path cluster");
                        }
                    }
                }
            }
            INFO(eulerian_clusters << " Eulerian clusters");
            INFO(clusters << " clusters");
            INFO(path_clusters << " path clusters.");
        }

        void AnalyzePathClusters(const string& stats_base_path, size_t distance) {
            INFO("Distance: " << distance);
            contracted_graph::ContractedGraphBuilder contracted_graph_builder(gp_.g, unique_storage_);
            INFO("Building contracted graph");

            ClusterGraphAnalyzerTester tester(gp_, contracted_graph_builder);
            tester.LaunchTests();

            auto contracted_graph = contracted_graph_builder.BuildContractedGraphFromDBG();
            auto scaffold_graph_constructor = scaffold_graph_utils::ScaffoldGraphConstructor(unique_storage_, distance, gp_.g);
            auto scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphUsingDijkstra();
            auto contracted_scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphFromContractedGraph(contracted_graph);
            INFO("Scaffold graph edges: " << scaffold_graph.EdgeCount());
            INFO("Contracted scaffold graph edges: " << contracted_scaffold_graph.EdgeCount());

            const size_t edge_read_threshold = configs_.edge_read_cluster_threshold;
            const size_t analyzer_read_threshold = configs_.global_read_cluster_threshold;
            auto cluster_storage_builder = cluster_storage::ClusterStorageBuilder(gp_.g, scaffold_graph,
                                                                                  barcode_extractor_ptr_, unique_storage_,
                                                                                  distance, edge_read_threshold);
            auto cluster_storage = cluster_storage_builder.ConstructClusterStorage();

            auto min_read_filter = make_shared<cluster_storage::MinReadClusterFilter>(analyzer_read_threshold);
            cluster_storage::ClusterStorageExtractor cluster_storage_extractor;
            auto high_covered_clusters = cluster_storage_extractor.FilterClusterStorage(cluster_storage, min_read_filter);

            cluster_storage::ClusterGraphAnalyzer ordering_analyzer(contracted_graph_builder);

            CountBasicClusterStats(cluster_storage, ordering_analyzer);

            INFO("Building reference transition storage")
            const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
            transitions::ContigPathBuilder contig_path_builder(gp_);
            auto transition_builder = make_shared<transitions::StrictTransitionStorageBuilder>(contig_path_builder);
            auto reference_transition_storage = BuildFilteredTransitionStorage(contig_path_builder, transition_builder,
                                                                               unique_storage_, path_to_reference);
            INFO("Reference transition storage size: " << reference_transition_storage.Size());

            INFO("Building path cluster storage")
            cluster_statistics::PathClusterStorageBuilder path_cluster_builder;
            auto path_cluster_storage =
                path_cluster_builder.BuildPathClusterStorage(ordering_analyzer, high_covered_clusters);
            INFO(path_cluster_storage.Size() << " distinct clusters");

            INFO("Analyzing path clusters")
            cluster_statistics::ClusterStorageAnalyzer cluster_analyzer(ordering_analyzer, scaffold_graph,
                                                                        reference_transition_storage,
                                                                        path_cluster_storage, high_covered_clusters);

            cluster_analyzer.FillStatistics();
            cluster_analyzer.SerializeStatistics(stats_base_path);

            INFO("Analyzing local contracted graph stats")
            contracted_graph::ContractedGraphClusterStatistics contracted_cluster_analyzer(gp_.g, contracted_graph,
                                                                                           ordering_analyzer,
                                                                                           high_covered_clusters,
                                                                                           reference_transition_storage);
            contracted_cluster_analyzer.FillStatistics();
            contracted_cluster_analyzer.SerializeStatistics(stats_base_path);
        }

        void AnalyzeContigPaths(const string& stats_base_path) {
            const string reference_path = cfg::get().ts_res.statistics.genome_path;
            const string cloud_contigs_path = cfg::get().ts_res.statistics.cloud_contigs_path;
            INFO("Reference path: " << reference_path);
            INFO("Cloud contig path: " << cloud_contigs_path);

            transitions::ContigPathBuilder contig_path_builder(gp_);
            transitions::StrictTransitionStorageBuilder strict_transition_builder(contig_path_builder);
            auto reference_paths = contig_path_builder.GetContigPaths(reference_path);
            auto contig_paths = contig_path_builder.GetContigPaths(cloud_contigs_path);

            const size_t min_length = 20000;
            path_extend::ScaffoldingUniqueEdgeStorage large_unique_storage;
            path_extend::ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, 100);
            unique_edge_analyzer.FillUniqueEdgeStorage(large_unique_storage);
            INFO(reference_paths.size() << " reference paths");
            INFO(contig_paths.size() << " contig paths");

            transitions::ContigPathFilter contig_path_filter(large_unique_storage);
            auto filtered_reference_paths = contig_path_filter.FilterPathsUsingUniqueStorage(reference_paths);
            auto filtered_contig_paths = contig_path_filter.FilterPathsUsingUniqueStorage(contig_paths);
            INFO(filtered_reference_paths.size() << " filtered reference paths");
            INFO(filtered_contig_paths.size() << " filtered contig paths");

            const string reference_output_path = fs::append_path(stats_base_path, "reference_paths");
            ofstream fout_reference(reference_output_path);
            PrintContigPaths(reference_paths, fout_reference);

            const string filtered_output_path = fs::append_path(stats_base_path, "filtered_reference_paths");
            ofstream fout_filt_ref(filtered_output_path);
            PrintContigPaths(filtered_reference_paths, fout_filt_ref);

            transitions::PathStatisticsExtractor path_stats_extractor(min_length, filtered_reference_paths,
                                                                      filtered_contig_paths, gp_.g);
            path_stats_extractor.FillStatistics();
            path_stats_extractor.SerializeStatistics(stats_base_path);

            INFO("Constructing scaffold graph");
            const size_t scaffolding_distance = 100000;
            auto scaffold_graph_constructor = scaffold_graph_utils::ScaffoldGraphConstructor(large_unique_storage,
                                                                                             scaffolding_distance, gp_.g);
            auto scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphUsingDijkstra();
            INFO(scaffold_graph.EdgeCount() << " edges in scaffold graph");
            const string scaffold_graph_path = fs::append_path(stats_base_path, "scaffold_graph");
            ofstream fout(scaffold_graph_path);
            scaffold_graph.Print(fout);


            INFO("Constructing initial tenx filter");
            auto tenx_resolver_configs = cfg::get().ts_res.tenx;
            typedef barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_t;
            auto barcode_extractor_ptr = make_shared<barcode_extractor_t>(gp_.barcode_mapper_ptr, gp_.g);
            path_extend::InitialTenXFilter initial_tenx_filter(gp_.g, barcode_extractor_ptr, tenx_resolver_configs);

            INFO("Constructing initial filter analyzer");
            transitions::InitialFilterStatisticsExtractor initial_filter_extractor(scaffold_graph,
                                                                                   reference_paths,
                                                                                   barcode_extractor_ptr,
                                                                                   initial_tenx_filter,
                                                                                   gp_.g, large_unique_storage);
            INFO("Analyzing initial filter");
            initial_filter_extractor.FillStatistics();
            initial_filter_extractor.SerializeStatistics(stats_base_path);
        }

        void PrintContigPaths(const vector<vector<transitions::EdgeWithMapping>>& contig_paths, ostream& stream) {
            size_t path_id = 0;
            for (const auto& path: contig_paths) {
                stream << "Path id: " << path_id << std::endl;
                stream << path.size() << " edges." << std::endl;
                for (const auto& edge: path) {
                    stream << edge.edge_.int_id() << std::endl;
                }
                ++path_id;
            }
        }

        void AnalyzeContractedGraph(const string& stats_base_path, size_t distance) {
            INFO("Distance: " << distance);
            contracted_graph::ContractedGraphBuilder contracted_graph_builder(gp_.g, unique_storage_);
            INFO("Building contracted graph");
            auto contracted_graph = contracted_graph_builder.BuildContractedGraphFromDBG();
            auto scaffold_graph_constructor = scaffold_graph_utils::ScaffoldGraphConstructor(unique_storage_, distance, gp_.g);
            auto scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphUsingDijkstra();
            auto contracted_scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphFromContractedGraph(contracted_graph);
            INFO("Scaffold graph edges: " << scaffold_graph.EdgeCount());
            INFO("Contracted scaffold graph edges: " << contracted_scaffold_graph.EdgeCount());
            const size_t edge_read_threshold = configs_.edge_read_cluster_threshold;
            const size_t analyzer_read_threshold = configs_.global_read_cluster_threshold;

            auto cluster_storage_builder = cluster_storage::ClusterStorageBuilder(gp_.g, scaffold_graph,
                                                                                  barcode_extractor_ptr_, unique_storage_,
                                                                                  distance, edge_read_threshold);
            auto cluster_storage = cluster_storage_builder.ConstructClusterStorage();


            cluster_storage::ClusterGraphAnalyzer ordering_analyzer(contracted_graph_builder);
            auto min_read_filter = make_shared<cluster_storage::MinReadClusterFilter>(analyzer_read_threshold);
            cluster_storage::ClusterStorageExtractor cluster_extractor;
            auto high_covered_clusters = cluster_extractor.FilterClusterStorage(cluster_storage, min_read_filter);

            cluster_statistics::PathClusterStorageBuilder path_cluster_builder;
            auto path_cluster_storage =
                path_cluster_builder.BuildPathClusterStorage(ordering_analyzer, high_covered_clusters);
            INFO(path_cluster_storage.Size() << " distinct clusters");

            auto name_to_transition_storage = BuildTransitionStorages(unique_storage_);
            auto reference_transition_storage = name_to_transition_storage.at("Reference");

            scaffold_graph_utils::ScaffoldGraphAnalyzer scaffold_analyzer(contracted_scaffold_graph);
            scaffold_analyzer.FillStatistics();
            scaffold_analyzer.SerializeStatistics(stats_base_path);

            contracted_graph::ContractedGraphAnalyzer contracted_analyzer(gp_.g, ordering_analyzer,
                                                                          *barcode_extractor_ptr_, path_cluster_storage,
                                                                          contracted_graph, name_to_transition_storage,
                                                                          reference_transition_storage,
                                                                          cluster_storage,
                                                                          analyzer_read_threshold);
            contracted_analyzer.FillStatistics();
            contracted_analyzer.SerializeStatistics(stats_base_path);
        }

        void LaunchAnalyzerForMultipleDistances(const string& stats_base_path,
                                                const std::function<void(const string&, size_t)> analyzer_function,
                                                const vector<size_t>& distances, size_t threads) {
#pragma omp parallel for num_threads(threads)
            for (size_t i = 0; i < distances.size(); ++i) {
                size_t distance = distances[i];
                string stat_path = fs::append_path(stats_base_path, "distance_" + std::to_string(distance));
                fs::make_dir(stat_path);
                INFO(stat_path);
                analyzer_function(stat_path, distance);
            }
        }
        DECL_LOGGER("StatisticsLauncher");
    };
}