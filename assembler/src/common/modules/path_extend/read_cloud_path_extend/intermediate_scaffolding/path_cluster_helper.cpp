//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_cluster_helper.hpp"

#include "simple_graph_utils.hpp"
#include "contracted_graph_from_simple.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/graph_cluster_storage_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/utils/barcode_score_functions.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

namespace path_extend {
namespace read_cloud {
std::vector<cluster_storage::Cluster> PathClusterExtractorHelper::GetPathClusters(
    const PathClusterExtractorHelper::SimpleTransitionGraph &graph) const {
    cluster_storage::GraphClusterStorageBuilder cluster_storage_builder(g_, barcode_extractor_, linkage_distance_);
    DEBUG("Constructing cluster storage");
    auto cluster_storage = cluster_storage_builder.ConstructClusterStorage(*initial_cluster_storage_, graph);
    ContractedGraphFromSimpleHelper contracted_helper(g_);
    cluster_storage::ClusterGraphAnalyzer cluster_graph_analyzer(contracted_helper);
    auto path_cluster_filter_ptr = std::make_shared<cluster_storage::PathClusterFilter>(cluster_graph_analyzer);
    cluster_storage::ClusterStorageExtractor cluster_extractor;
    DEBUG("Processing");
    auto path_clusters = cluster_extractor.FilterClusterStorage(cluster_storage, path_cluster_filter_ptr);

    return path_clusters;
}
PathClusterExtractorHelper::PathClusterExtractorHelper(
        const Graph &g,
        std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t linkage_distance)
    : g_(g),
      initial_cluster_storage_(initial_cluster_storage),
      barcode_extractor_(barcode_extractor),
      linkage_distance_(linkage_distance) {}

PathClusterStorage GraphBasedPathClusterNormalizer::GetNormalizedStorage(
    const std::vector<cluster_storage::Cluster> &path_clusters) const {
    typedef std::set<scaffold_graph::ScaffoldVertex> VertexSet;
    std::map<VertexSet, double> cluster_to_weight;

    for (const auto &cluster: path_clusters) {
        const auto vertex_set = cluster.GetVertexSet();
        if (cluster_to_weight.find(vertex_set) == cluster_to_weight.end()) {
            cluster_to_weight[vertex_set] = 1;
        } else {
            cluster_to_weight[vertex_set]++;
        }
    }

    size_t local_tail_threshold = 10000;
    for (auto &entry: cluster_to_weight) {
        if (entry.first.size() == 2) {
            double length_normalizer = 0;
            for (const auto &vertex: entry.first) {
                length_normalizer += static_cast<double>(std::max(local_tail_threshold, vertex.GetLengthFromGraph(g_)));
            }
            entry.second = entry.second / length_normalizer * 100000;
            for (const auto &vertex: entry.first) {
                entry.second /= vertex.GetCoverageFromGraph(g_);
            }
            entry.second *= 10000;
        }
    }
    PathClusterStorage storage(cluster_to_weight);
    return storage;
}
GraphBasedPathClusterNormalizer::GraphBasedPathClusterNormalizer(const Graph &g_) : g_(g_) {}

PathClusterStorage::PathClusterStorage(const std::map<PathClusterStorage::VertexSet, double> &cluster_to_weight)
    : cluster_to_weight_(cluster_to_weight) {}
PathClusterStorage::ClusterToWeightT::const_iterator PathClusterStorage::begin() const {
    return cluster_to_weight_.begin();
}
PathClusterStorage::ClusterToWeightT::const_iterator PathClusterStorage::end() const {
    return cluster_to_weight_.end();
}

bool PathClusterConflictResolver::AreClustersConflicted(
    const PathClusterConflictResolver::VertexSet &first,
    const PathClusterConflictResolver::VertexSet &second,
    const PathClusterConflictResolver::SimpleTransitionGraph &graph) const {
    TRACE("First size: " << first.size());
    TRACE("Second size: " << second.size());
    //todo maybe it is not needed
    if (first.size() != second.size() or first.size() == 1) {
        return false;
    }

    VertexSet intersection;
    std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                          std::inserter(intersection, intersection.begin()));
    size_t possible_overlap = intersection.size();
    if (possible_overlap != first.size() - 1) {
        return false;
    }
    cluster_storage::GraphAnalyzer graph_analyzer;

    SubgraphGetter subgraph_getter;
    auto first_subgraph = subgraph_getter.GetSubgraph(graph, first);
    auto second_subgraph = subgraph_getter.GetSubgraph(graph, second);
    std::vector<std::vector<ScaffoldVertex>> first_ham_paths = graph_analyzer.GetHamiltonianPaths(first_subgraph);
    std::vector<std::vector<ScaffoldVertex>> second_ham_paths = graph_analyzer.GetHamiltonianPaths(second_subgraph);

    for (const auto &first_path: first_ham_paths) {
        for (const auto &second_path: second_ham_paths) {
            if (CheckOverlap(first_path, second_path, possible_overlap) or
                CheckOverlap(second_path, first_path, possible_overlap)) {
                return false;
            }
        }
    }
    return true;
}

std::vector<PathClusterConflictResolver::VertexSet> PathClusterConflictResolver::GetClusterSets(
    const PathClusterConflictResolver::SimpleTransitionGraph &graph, const PathClusterStorage &storage) const {
    std::set<VertexSet> result_set;
    //fixme move from here
    const size_t MAX_CLOUD_SIZE = 6;
    DEBUG("Checking entries");
    for (const auto &entry: storage) {
        if (entry.first.size() <= MAX_CLOUD_SIZE) {
            result_set.insert(entry.first);
            for (const auto &vertex: entry.first) {
                VERIFY_DEV(graph.ContainsVertex(vertex));
            }
        }
    }
    std::set<ScaffoldVertex> vertices;
    for (const auto &vertex: graph) {
        vertices.insert(vertex);
        vertices.insert(vertex.GetConjugateFromGraph(g_));
    }
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = 10000;
    size_t count_threshold = 1;
    const size_t max_threads = 1;
    auto scaffold_vertex_index = helper.TailEdgeScaffoldVertexIndex(g_, *barcode_extractor_, vertices,
                                                                    count_threshold, tail_threshold, max_threads);
    VERIFY_DEV(math::gr(relative_threshold_, 0.));
    size_t path_cluster_size = storage.cluster_to_weight_.size();
    DEBUG(path_cluster_size << " path clusters");
    size_t processed = 0;
    for (const auto &entry: storage) {
        const auto &first_set = entry.first;
        for (const auto &other: storage) {
            const auto &second_set = other.first;
            TRACE("Checking conflicts");
            if (AreClustersConflicted(first_set, second_set, graph)) {
                DEBUG("Found conflict");
                double clash_score = GetClashScore(first_set, second_set, scaffold_vertex_index);
                DEBUG("Got clash score");
                if (math::ge(clash_score, relative_threshold_)) {
                    result_set.erase(second_set);
                } else if (math::le(clash_score, 1 / relative_threshold_)) {
                    result_set.erase(first_set);
                } else {
                    result_set.erase(first_set);
                    result_set.erase(second_set);
                }
            }
            ++processed;
            if (path_cluster_size > 10 and processed % (path_cluster_size * path_cluster_size / 10) == 0) {
                DEBUG("Processed " << processed << " pairs");
            }
        }
    }
    DEBUG("Returning");
    std::vector<VertexSet> result;
    std::move(result_set.begin(), result_set.end(), std::back_inserter(result));
    return result;
}
PathClusterConflictResolver::PathClusterConflictResolver(
        const Graph &g,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        double relative_threshold) :
    g_(g), barcode_extractor_(barcode_extractor), relative_threshold_(relative_threshold) {}
bool PathClusterConflictResolver::CheckOverlap(const std::vector<PathClusterConflictResolver::ScaffoldVertex> &first,
                                               const std::vector<PathClusterConflictResolver::ScaffoldVertex> &second,
                                               size_t overlap) const {
    VERIFY_DEV(first.size() >= overlap);
    for (size_t i = 0; i < overlap; ++i) {
        if (first[first.size() - overlap + i] != second[i]) {
            return false;
        }
    }
    return true;
}
double PathClusterConflictResolver::GetClashScore(const std::set<PathClusterConflictResolver::ScaffoldVertex> &first,
                                                  const std::set<PathClusterConflictResolver::ScaffoldVertex> &second,
                                                  PathClusterConflictResolver::ScaffoldBarcodeIndex scaffold_vertex_index) const {
    //todo make one pass instead of three
    //todo deal with conjugates
    std::set<ScaffoldVertex> intersection;
    std::set<ScaffoldVertex> first_minus_second;
    std::set<ScaffoldVertex> second_minus_first;
    std::set_intersection(first.begin(), first.end(), second.begin(), second.end(),
                          std::inserter(intersection, intersection.end()));
    std::set_difference(first.begin(), first.end(), second.begin(), second.end(),
                        std::inserter(first_minus_second, first_minus_second.end()));
    std::set_difference(first.begin(), first.end(), second.begin(), second.end(),
                        std::inserter(second_minus_first, second_minus_first.end()));
    std::set<ScaffoldVertex> first_minus_second_inter;
    std::set<ScaffoldVertex> second_minus_first_inter;
    auto intersection_barcodes = GetBarcodesFromSet(intersection, scaffold_vertex_index);
    auto first_minus_second_barcodes = GetBarcodesFromSet(first_minus_second, scaffold_vertex_index);
    auto second_minus_first_barcodes = GetBarcodesFromSet(second_minus_first, scaffold_vertex_index);
    ContainmentIndex containment_index;
    TRACE("Getting scores");
    double
        first_minus_second_ci = containment_index.GetScoreFromSets(first_minus_second_barcodes, intersection_barcodes);
    double
        second_minus_first_ci = containment_index.GetScoreFromSets(second_minus_first_barcodes, intersection_barcodes);
    TRACE("Got scores");
    const double infty = relative_threshold_ + 1;
    if (math::eq(second_minus_first_ci, 0.)) {
        if (math::eq(first_minus_second_ci, 0.)) {
            return 1;
        } else {
            return infty;
        }
    }
    return first_minus_second_ci / second_minus_first_ci;
}
std::set<barcode_index::BarcodeId> PathClusterConflictResolver::GetBarcodesFromSet(
    const std::set<PathClusterConflictResolver::ScaffoldVertex> &vertices,
    PathClusterConflictResolver::ScaffoldBarcodeIndex scaffold_vertex_index) const {
    std::set<barcode_index::BarcodeId> result;
    for (const auto &vertex: vertices) {
        auto head_barcodes_begin = scaffold_vertex_index->GetHeadBegin(vertex);
        auto head_barcodes_end = scaffold_vertex_index->GetHeadEnd(vertex);
        auto tail_barcodes_begin = scaffold_vertex_index->GetTailBegin(vertex);
        auto tail_barcodes_end = scaffold_vertex_index->GetTailEnd(vertex);
        std::copy(head_barcodes_begin, head_barcodes_end, std::inserter(result, result.end()));
        std::copy(tail_barcodes_begin, tail_barcodes_end, std::inserter(result, result.end()));
    }
    return result;
}

std::vector<CloudPathExtractor::InternalPathWithSet> CloudPathExtractor::ExtractAllPaths(
    const CloudPathExtractor::SimpleTransitionGraph &graph,
    const ScaffoldVertex &source, const ScaffoldVertex &sink) const {
    DEBUG("Extracting all paths");
    TRACE(source.int_id() << " -> " << sink.int_id());
    std::vector<InternalPathWithSet> result;
    std::queue<InternalPathWithSet> current_paths;
    InternalPathWithSet start;
    start.AddVertex(source);
    current_paths.push(start);
    size_t ops_counter = 0;
    while (not current_paths.empty()) {
        auto last_path = current_paths.front();
        string last_path_string;
        TRACE("Last path size: " << last_path.path_.size());
        for (const auto &vertex: last_path.path_) {
            last_path_string += std::to_string(vertex.int_id()) += " -> ";
            TRACE("Last path:");
            TRACE(last_path_string);
        }
        current_paths.pop();
        ScaffoldVertex last_vertex = last_path.path_.back();
        TRACE("Last vertex in graph: " << graph.ContainsVertex(last_vertex));
        if (last_vertex == sink) {
            result.push_back(last_path);
            continue;
        }
        for (const auto &new_vertex: graph.OutNeighbours(last_vertex)) {
            TRACE("New vertex: " << new_vertex.int_id());
            if (not last_path.HasVertex(new_vertex)) {
                auto last_path_copy = last_path;
                last_path_copy.AddVertex(new_vertex);
                TRACE("Adding new path");
                current_paths.push(last_path_copy);
            }
        }
        ++ops_counter;
        if (ops_counter > 200) {
            DEBUG("Too many operations!");
            std::vector<InternalPathWithSet> empty;
            return empty;
        }
    }
    DEBUG("Extracted " << result.size() << " paths");
    return result;
}
bool CloudPathExtractor::IsPathCorrect(const InternalPathWithSet &path,
                                       const std::vector<CloudPathExtractor::VertexSet> &clouds) const {
    DEBUG("Checking path");
    const auto vertex_path = path.path_;
    const auto &vertex_set = path.path_vertices_;
    std::vector<std::vector<ScaffoldVertex>> supporting_clouds;
    std::unordered_map<ScaffoldVertex, size_t> vertex_to_index;
    for (size_t i = 0; i < vertex_path.size(); ++i) {
        vertex_to_index[vertex_path[i]] = i;
    }

    auto path_compare_strong = [&vertex_to_index](const ScaffoldVertex &first, const ScaffoldVertex &second) {
      return vertex_to_index.at(first) < vertex_to_index.at(second);
    };

    auto path_compare_weak = [&vertex_to_index](const ScaffoldVertex &first, const ScaffoldVertex &second) {
      return vertex_to_index.at(first) <= vertex_to_index.at(second);
    };

    for (const auto &cloud: clouds) {
        bool cloud_in_path = true;
        for (const auto &vertex: cloud) {
            if (vertex_set.find(vertex) == vertex_set.end()) {
                cloud_in_path = false;
                break;
            }

        }
        if (not cloud_in_path) {
            continue;
        }
        std::vector<ScaffoldVertex> cloud_vertices;
        std::copy(cloud.begin(), cloud.end(), std::back_inserter(cloud_vertices));
        std::sort(cloud_vertices.begin(), cloud_vertices.end(), path_compare_strong);
        for (size_t i = 1; i < cloud_vertices.size(); ++i) {
            if (vertex_to_index.at(cloud_vertices[i]) != vertex_to_index.at(cloud_vertices[i - 1]) + 1) {
                break;
            }
        }
        supporting_clouds.push_back(cloud_vertices);
    }
    DEBUG("Found " << supporting_clouds.size() << " supporting clouds");
    const size_t max_supporting = 20;
    if (supporting_clouds.size() > max_supporting) {
        return false;
    }

    for (size_t i = 0; i < vertex_path.size(); ++i) {
        for (size_t j = i + 1; j < vertex_path.size(); ++j) {
            if (i != 0 or j != vertex_path.size() - 1) {
                bool found_intersection = false;
                DEBUG("Checking intersections");
                for (const auto &cloud: supporting_clouds) {
                    ScaffoldVertex cloud_first = cloud[0];
                    ScaffoldVertex cloud_last = cloud.back();
                    ScaffoldVertex path_first = vertex_path[i];
                    ScaffoldVertex path_last = vertex_path[j];
                    bool right_nontrival = path_compare_weak(cloud_first, path_last) and
                        path_compare_strong(path_first, cloud_first) and path_compare_strong(path_last, cloud_last);
                    bool left_nontrivial = path_compare_weak(path_first, cloud_last) and
                        path_compare_strong(cloud_last, path_last) and path_compare_strong(cloud_first, path_first);
                    if (left_nontrivial or right_nontrival) {
                        found_intersection = true;
                    }
                }
                DEBUG("Checked intersections");
                if (not(found_intersection)) {
                    return false;
                }
            } else if (j - i <= 2) {
                bool found_containing = false;
                for (const auto &cloud: supporting_clouds) {
                    if (cloud[0] == vertex_path[i] and cloud.back() == vertex_path[j]) {
                        found_containing = true;
                    }
                }
                if (not(found_containing)) {
                    return false;
                }
            }
        }
    }
    return true;
}
std::vector<std::vector<CloudPathExtractor::ScaffoldVertex>> CloudPathExtractor::ExtractCorrectPaths(
    const CloudPathExtractor::SimpleTransitionGraph &graph,
    const CloudPathExtractor::ScaffoldVertex &source,
    const CloudPathExtractor::ScaffoldVertex &sink,
    const std::vector<CloudPathExtractor::VertexSet> &clouds) const {
    DEBUG("Extracting correct paths");
    const auto paths = ExtractAllPaths(graph, source, sink);
    DEBUG("Extracted all paths");
    std::vector<std::vector<ScaffoldVertex>> result;
    if (paths.size() == 1) {
        DEBUG("Single path");
        result.push_back(paths[0].path_);
        return result;
    }
    if (paths.size() > 20) {
        DEBUG("Too many paths");
        return result;
    }
    for (const auto &path: paths) {
        if (IsPathCorrect(path, clouds)) {
            result.push_back(path.path_);
        }
    }
    return result;
}
void CloudPathExtractor::InternalPathWithSet::AddVertex(const CloudPathExtractor::ScaffoldVertex &vertex) {
    VERIFY_DEV(not HasVertex(vertex));
    path_.push_back(vertex);
    path_vertices_.insert(vertex);
}
bool CloudPathExtractor::InternalPathWithSet::HasVertex(const CloudPathExtractor::ScaffoldVertex &vertex) const {
    return path_vertices_.find(vertex) != path_vertices_.end();
}
std::vector<std::vector<CloudBasedPathExtractor::ScaffoldVertex>> CloudBasedPathExtractor::GetCorrectPaths(
    const CloudBasedPathExtractor::SimpleTransitionGraph &graph,
    const CloudBasedPathExtractor::ScaffoldVertex &source,
    const CloudBasedPathExtractor::ScaffoldVertex &sink) const {
    PathClusterExtractorHelper path_cluster_extractor(g_, initial_cluster_storage_,
                                                      barcode_extractor_, linkage_distance_);
    DEBUG("Getting path clusters");
    auto path_clusters = path_cluster_extractor.GetPathClusters(graph);
    GraphBasedPathClusterNormalizer path_cluster_normalizer(g_);
    DEBUG("Get normalized storage")
    auto cluster_to_weight = path_cluster_normalizer.GetNormalizedStorage(path_clusters);

    PathClusterConflictResolver conflict_resolver(g_, barcode_extractor_, relative_cluster_threshold_);
    DEBUG("Resolving conflicts")
    auto final_clusters = conflict_resolver.GetClusterSets(graph, cluster_to_weight);
    CloudPathExtractor cloud_path_extractor;
    DEBUG("Extracting paths using final clusters");
    auto resulting_paths = cloud_path_extractor.ExtractCorrectPaths(graph, source, sink, final_clusters);
    DEBUG("Extracted paths");
    return resulting_paths;
}
CloudBasedPathExtractor::CloudBasedPathExtractor(
    const Graph &g,
    std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t linkage_distance, double relative_cluster_threshold)
    : g_(g),
      initial_cluster_storage_(initial_cluster_storage),
      barcode_extractor_(barcode_extractor),
      linkage_distance_(linkage_distance),
      relative_cluster_threshold_(relative_cluster_threshold) {}

std::vector<cluster_storage::Cluster> ScaffoldGraphPathClusterHelper::GetPathClusters(
    const scaffold_graph::ScaffoldGraph &graph) const {
    const size_t linkage_distance = 10;
    TransitionGraph simple_graph;
    for (const auto &vertex: graph.vertices()) {
        simple_graph.AddVertex(vertex);
    }
    for (const auto &edge: graph.edges()) {
        simple_graph.AddEdge(edge.getStart(), edge.getEnd());
    }
    PathClusterExtractorHelper cluster_extractor_helper(g_, initial_cluster_storage_,
                                                        barcode_extractor_, linkage_distance);
    return cluster_extractor_helper.GetPathClusters(simple_graph);
}
ScaffoldGraphPathClusterHelper::ScaffoldGraphPathClusterHelper(
    const Graph &g,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
    size_t max_threads)
    : g_(g), barcode_extractor_(barcode_extractor),
      initial_cluster_storage_(initial_cluster_storage), max_threads_(max_threads) {}

std::vector<cluster_storage::Cluster> ScaffoldGraphPathClusterHelper::GetAllClusters(
    const scaffold_graph::ScaffoldGraph &graph) const {
    const size_t linkage_distance = 10;
    cluster_storage::GraphClusterStorageBuilder cluster_storage_builder(g_, barcode_extractor_, linkage_distance);
    DEBUG("Constructing cluster storage");
    TransitionGraph simple_graph;
    for (const auto &vertex: graph.vertices()) {
        simple_graph.AddVertex(vertex);
    }
    for (const auto &edge: graph.edges()) {
        simple_graph.AddEdge(edge.getStart(), edge.getEnd());

    }
    auto cluster_storage = cluster_storage_builder.ConstructClusterStorage(*initial_cluster_storage_, simple_graph);
    std::vector<cluster_storage::Cluster> result;
    for (const auto &entry: cluster_storage) {
        if (entry.second.Size() >= 2) {
            result.push_back(entry.second);
        }
    }
    return result;
}
std::vector<ScaffoldGraphPathClusterHelper::Cluster> ScaffoldGraphPathClusterHelper::GetPathClusters(
    const std::vector<ScaffoldGraphPathClusterHelper::Cluster> &clusters) const {
    ContractedGraphFromSimpleHelper contracted_helper(g_);
    cluster_storage::ClusterGraphAnalyzer cluster_graph_analyzer(contracted_helper);
    auto path_cluster_filter_ptr = std::make_shared<cluster_storage::PathClusterFilter>(cluster_graph_analyzer);
    cluster_storage::ClusterStorageExtractor cluster_extractor;
    return cluster_extractor.FilterClusters(clusters, path_cluster_filter_ptr);
}
std::vector<std::set<ScaffoldGraphPathClusterHelper::ScaffoldVertex>> ScaffoldGraphPathClusterHelper::GetCorrectedClusters(
    const std::vector<ScaffoldGraphPathClusterHelper::Cluster> &path_clusters,
    const scaffold_graph::ScaffoldGraph &graph) const {
    auto transition_graph = ScaffoldToTransition(graph);
    return GetCorrectedClusters(path_clusters, transition_graph);
}
std::vector<std::set<ScaffoldGraphPathClusterHelper::ScaffoldVertex>> ScaffoldGraphPathClusterHelper::GetFinalClusters(
    const scaffold_graph::ScaffoldGraph &graph) const {
    auto transition_graph = ScaffoldToTransition(graph);
    return GetFinalClusters(transition_graph);
}
ScaffoldGraphPathClusterHelper::TransitionGraph ScaffoldGraphPathClusterHelper::ScaffoldToTransition(
    const scaffold_graph::ScaffoldGraph &graph) const {
    TransitionGraph simple_graph;
    for (const auto &vertex: graph.vertices()) {
        simple_graph.AddVertex(vertex);
    }
    for (const auto &edge: graph.edges()) {
        simple_graph.AddEdge(edge.getStart(), edge.getEnd());
    }
    return simple_graph;
}
std::vector<std::set<ScaffoldGraphPathClusterHelper::ScaffoldVertex>> ScaffoldGraphPathClusterHelper::GetFinalClusters(
    const ScaffoldGraphPathClusterHelper::TransitionGraph &graph) const {
    const size_t linkage_distance = 10;
    PathClusterExtractorHelper cluster_extractor_helper(g_, initial_cluster_storage_,
                                                        barcode_extractor_, linkage_distance);
    auto path_clusters = cluster_extractor_helper.GetPathClusters(graph);
    return GetCorrectedClusters(path_clusters, graph);
}
std::vector<std::set<ScaffoldGraphPathClusterHelper::ScaffoldVertex>> ScaffoldGraphPathClusterHelper::GetCorrectedClusters(
    const std::vector<ScaffoldGraphPathClusterHelper::Cluster> &path_clusters,
    const ScaffoldGraphPathClusterHelper::TransitionGraph &graph) const {
    GraphBasedPathClusterNormalizer path_cluster_normalizer(g_);
    auto cluster_to_weight = path_cluster_normalizer.GetNormalizedStorage(path_clusters);
    const double relative_threshold = 2;
    PathClusterConflictResolver conflict_resolver(g_, barcode_extractor_, relative_threshold);
    auto final_clusters = conflict_resolver.GetClusterSets(graph, cluster_to_weight);
    return final_clusters;
}
}
}