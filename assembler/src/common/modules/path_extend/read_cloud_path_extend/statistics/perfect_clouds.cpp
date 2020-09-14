//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_clouds.hpp"

#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_helper.hpp"

namespace path_extend {
namespace read_cloud {

PerfectClustersAnalyzer::SetDistribution PerfectClustersAnalyzer::ConstructPerfectClusters(
    const scaffold_graph::ScaffoldGraph &perfect_graph) const {

    const size_t min_read_threshold = 1;
    const auto &g = gp_.get<Graph>();
    const auto &barcode_mapper = gp_.get<barcode_index::FrameBarcodeIndex<Graph>>();
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper, g);
    cluster_storage::HalfEdgeClusterStorageHelper cluster_storage_helper(g, barcode_extractor,
                                                                         min_read_threshold, max_threads_);
    std::set<scaffold_graph::ScaffoldVertex> vertices;
    std::copy(perfect_graph.vbegin(), perfect_graph.vend(), std::inserter(vertices, vertices.end()));
    auto initial_storage_builder = cluster_storage_helper.GetInitialStorageBuilder(vertices);
    DEBUG("Constructing initial storage");
    auto initial_storage = std::make_shared<cluster_storage::InitialClusterStorage>(
        initial_storage_builder->ConstructInitialClusterStorage());

    ScaffoldGraphPathClusterHelper path_cluster_helper(g, barcode_extractor, initial_storage, max_threads_);
    std::vector<Cluster> raw_clusters = path_cluster_helper.GetAllClusters(perfect_graph);
    SetDistribution result;
    for (const auto &cluster: raw_clusters) {
        VertexSet cluster_vertices;
        for (const auto &entry: cluster) {
            cluster_vertices.insert(entry.first);
        }
        result[cluster_vertices] += 1;
    }
    return result;
}
void PerfectClustersAnalyzer::AnalyzePerfectClouds(const std::string &path_to_reference, size_t min_length) const {
    const auto &g = gp_.get<Graph>();
    PerfectScaffoldGraphConstructor constructor(g);
    omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(g);
    std::vector<EdgeId> long_edges;
    for (const auto &edge: edge_it_helper) {
        if (g.length(edge) >= min_length) {
            long_edges.push_back(edge);
        }
    }

    std::vector<size_t> length_thresholds;
    const size_t max_length = 20000;
    const size_t step = 1000;
    for (size_t length = min_length; length <= max_length; length += step) {
        length_thresholds.push_back(length);
    }
    validation::FilteredReferencePathHelper path_helper(gp_.get<Graph>(),
                                                        gp_.get<EdgeIndex<Graph>>(),
                                                        gp_.get<KmerMapper<Graph>>());
    ScaffoldingUniqueEdgeAnalyzer unique_analyzer(gp_, min_length, 5000.0);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_analyzer.FillUniqueEdgeStorage(unique_storage);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromUnique(path_to_reference, unique_storage);
    auto perfect_scaffold_graph = constructor.ConstuctPerfectGraph(reference_paths.paths_, min_length);
    auto perfect_clusters = ConstructPerfectClusters(perfect_scaffold_graph);

    std::string output_path = fs::append_path(output_dir_, "perfect_cluster_stats");
    std::ofstream fout(output_path);
    fout << "length total_edges mean_in_cluster\n";
    for (const auto &length_threshold: length_thresholds) {

        INFO("Length threshold: " << length_threshold);
        size_t total_edges = 0;
        for (const auto &edge: long_edges) {
            if (g.length(edge) >= length_threshold) {
                ++total_edges;
            }
        }
        INFO("Total edges: " << total_edges);
        double mean_edges_num = GetMeanEdgeNumber(perfect_clusters, length_threshold, g);
        INFO("Mean edges: " << mean_edges_num);
        fout << length_threshold << " " << total_edges << " " << mean_edges_num << "\n";
    }
}
double PerfectClustersAnalyzer::GetMeanEdgeNumber(const SetDistribution &clusters,
                                                  size_t length_threshold, const Graph &g) const {
    size_t total_edges = 0;
    size_t total_clusters = 0;
    for (const auto &entry: clusters) {
        const auto &vertex_set = entry.first;
        size_t long_edges = std::count_if(vertex_set.begin(), vertex_set.end(),
                                          [&g, length_threshold](const scaffold_graph::ScaffoldVertex &vertex) {
                                            return vertex.GetLengthFromGraph(g) >= length_threshold;
                                          });
        if (long_edges > 0) {
            total_edges += long_edges * entry.second;
            total_clusters += entry.second;
        }
    }
    if (total_clusters == 0) {
        return 0;
    }
    INFO("Total edges: " << total_edges);
    INFO("Total clusters: " << total_clusters);
    return static_cast<double>(total_edges) / static_cast<double>(total_clusters);
}
PerfectClustersAnalyzer::PerfectClustersAnalyzer(const GraphPack &gp,
                                                 const std::string &output_dir,
                                                 size_t max_threads) :
    gp_(gp),
    output_dir_(output_dir),
    max_threads_(max_threads) {}
}
}