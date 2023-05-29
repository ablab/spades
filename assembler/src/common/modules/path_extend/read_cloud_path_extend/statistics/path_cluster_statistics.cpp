#include "path_cluster_statistics.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "common/barcode_index/cluster_storage/initial_cluster_storage_builder.hpp"

namespace path_extend {

vector<SubgraphInfo> PathClusterStatisticsExtractor::GetAllSubgraphInfo(const ScaffoldGraphStorage &storage) {
    VERIFY_DEV(cfg::get().ts_res.debug_mode);
    vector<SubgraphInfo> result;
    ScaffoldGraphGapCloserParamsConstructor params_constructor;
    auto subgraph_extractor_params =
        params_constructor.ConstructSubgraphExtractorParamsFromConfig(storage.GetLargeLengthThreshold());
    auto path_extractor_params = params_constructor.ConstructPathExtractorParamsFromConfig();
    ScaffoldIndexInfoExtractorHelper scaffold_index_helper;
    auto scaffold_index_extractor =
        scaffold_index_helper.ConstructIndexExtractorFromParams(storage.GetSmallScaffoldGraph(), gp_,
                                                                subgraph_extractor_params);

    const size_t linkage_distance = path_extractor_params.linkage_distance_;
    const size_t min_read_threshold = path_extractor_params.min_read_threshold_;
    std::set<path_extend::scaffold_graph::ScaffoldVertex> target_edges;
    std::copy(storage.GetSmallScaffoldGraph().vbegin(), storage.GetSmallScaffoldGraph().vend(),
              std::inserter(target_edges, target_edges.begin()));
    DEBUG(target_edges.size() << " target edges.");
    auto barcode_extractor_ptr =
        make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    size_t cluster_storage_builder_threads = cfg::get().max_threads;
    auto edge_cluster_extractor =
        make_shared<cluster_storage::AccurateEdgeClusterExtractor>(gp_.g, barcode_extractor_ptr,
                                                                   linkage_distance, min_read_threshold);
    auto cluster_storage_builder =
        std::make_shared<cluster_storage::EdgeInitialClusterStorageBuilder>(gp_.g, edge_cluster_extractor,
                                                                            target_edges, linkage_distance,
                                                                            min_read_threshold,
                                                                            cluster_storage_builder_threads);
    INFO("Constructing initial cluster storage");
    auto initial_cluster_storage = std::make_shared<cluster_storage::InitialClusterStorage>(
        cluster_storage_builder->ConstructInitialClusterStorage());
    INFO("Constructed initial cluster storage");
    ScaffoldGraphExtractor scaffold_graph_extractor;
    auto univocal_edges = scaffold_graph_extractor.ExtractUnivocalEdges(storage.GetLargeScaffoldGraph());
    CloudScaffoldSubgraphExtractor subgraph_extractor(gp_.g, scaffold_index_extractor, subgraph_extractor_params);
    PathClusterExtractionParams path_cluster_extraction_params {gp_.g, initial_cluster_storage,
                                                                barcode_extractor_ptr, linkage_distance};
    validation::SimpleTransitionGraphValidatorConstructor validator_constructor(gp_);
    auto validator = validator_constructor.GetValidator(cfg::get().ts_res.statistics.genome_path);
    const bool reference_validation_on = true;

    for (const auto& edge: univocal_edges) {
        auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(storage.GetSmallScaffoldGraph(),
                                                                          edge.getStart(), edge.getEnd());
        auto pair_to_length = ConstructLengthMap(subgraph, storage.GetSmallScaffoldGraph());
        if (CheckSubgraph(subgraph, edge.getStart(), edge.getEnd(), validator, reference_validation_on)) {
            auto info =
                GetSubgraphInfo(subgraph, edge.getStart(), edge.getEnd(), pair_to_length, validator,
                                path_cluster_extraction_params, reference_validation_on);
            result.push_back(info);
        }
    }
    INFO("Analyzed " << result.size() << " subgraphs");
    return result;
}

SubgraphInfo PathClusterStatisticsExtractor::GetSubgraphInfo(
        const PathClusterStatisticsExtractor::SimpleTransitionGraph &graph,
        const scaffold_graph::ScaffoldVertex &source, const scaffold_graph::ScaffoldVertex &sink,
        const PathClusterStatisticsExtractor::ScaffoldEdgeMap &scaffold_edge_to_len,
        const validation::SimpleTransitionGraphValidator &validator,
        const PathClusterExtractionParams &path_cluster_extraction_params,
        bool reference_validation_on) const {
    PathClusterExtractorHelper path_cluster_extractor_helper(path_cluster_extraction_params.g_,
                                                             path_cluster_extraction_params.init_cluster_storage_,
                                                             path_cluster_extraction_params.barcode_extractor_,
                                                             path_cluster_extraction_params.linkage_distance_);
    DEBUG("Source: " << source.int_id());
    DEBUG("Sink: " << sink.int_id());
    DEBUG("Printing graph");
    for (const auto& vertex: graph) {
        for(auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            TRACE(vertex.int_id() << " -> " << (*it).int_id());
        }
    }
    DEBUG("Printed graph");

    auto path_clusters = path_cluster_extractor_helper.GetPathClusters(graph);
    INFO(path_clusters.size() << " path clusters");

    GraphBasedPathClusterNormalizer path_cluster_normalizer(gp_.g);
    auto cluster_to_weight = path_cluster_normalizer.GetNormalizedStorage(path_clusters);
    INFO("Normalized");
    const double relative_threshold = 2;
    PathClusterConflictResolver conflict_resolver(relative_threshold);
    auto final_clusters = conflict_resolver.GetClusterSets(graph, cluster_to_weight);
    INFO("Resolved conflicts");

    CloudPathExtractor correct_path_extractor;
    auto all_path_sets = correct_path_extractor.ExtractAllPaths(graph, source, sink);
    vector<vector<ScaffoldVertex>> all_paths;
    for (const auto& path_set: all_path_sets) {
        all_paths.push_back(path_set.path_);
    }
    auto resulting_paths = correct_path_extractor.ExtractCorrectPaths(graph, source, sink, final_clusters);

    vector<ScaffoldVertex> correct_path;
    if (reference_validation_on) {
        const auto correct_path_result = validator.GetCorrectPath(graph, source, sink);
        VERIFY_DEV(correct_path_result.is_initialized());
        correct_path = correct_path_result.get();
    }
    auto id_map = GetIdMap(graph, correct_path);
    std::map<ScaffoldVertex, size_t> vertex_to_len;
    std::map<ScaffoldVertex, double> vertex_to_cov;
    for (const auto& vertex: graph) {
        vertex_to_len[vertex] = vertex.getLengthFromGraph(gp_.g);
        vertex_to_cov[vertex] = vertex.getCoverageFromGraph(gp_.g);
    }
    SubgraphInfo result(graph, source, sink, cluster_to_weight, final_clusters, resulting_paths, all_paths,
                        correct_path, id_map, vertex_to_cov, vertex_to_len, scaffold_edge_to_len);
    return result;
}

PathClusterStatisticsExtractor::PathClusterStatisticsExtractor(const conj_graph_pack &gp) : gp_(gp) {}
string PathClusterStatisticsExtractor::RequestCorrectPath(const PathClusterStatisticsExtractor::SimpleTransitionGraph &graph,
                                                          const scaffold_graph::ScaffoldVertex &source,
                                                          const scaffold_graph::ScaffoldVertex &sink,
                                                          const validation::SimpleTransitionGraphValidator &validator) const {
    string path_string = "";
    boost::optional<vector<scaffold_graph::ScaffoldVertex>> correct_path =
        validator.GetCorrectPath(graph, source, sink);
    if (not correct_path.is_initialized()) {
        path_string = "No correct path!";
    } else {
        for (const auto& vertex: correct_path.get()) {
            path_string += std::to_string(vertex.int_id()) + " -> ";
        }
    }
    return path_string;
}
bool PathClusterStatisticsExtractor::CheckSubgraph(const PathClusterStatisticsExtractor::SimpleTransitionGraph &graph,
                                                   const scaffold_graph::ScaffoldVertex &source,
                                                   const scaffold_graph::ScaffoldVertex &sink,
                                                   const validation::SimpleTransitionGraphValidator &validator,
                                                   bool reference_validation_on) const {
    const size_t subgraph_size_threshold = 4;
    bool is_graph_large = graph.size() >= subgraph_size_threshold and graph.GetEdgesCount() >= graph.size();
    bool has_source_and_sink = graph.ContainsVertex(source) and graph.ContainsVertex(sink);
    if (not (is_graph_large and has_source_and_sink)) {
//        INFO("Subgraph is too small");
        return false;
    }

    if (reference_validation_on) {
        auto correct_path_result = validator.GetCorrectPath(graph, source, sink);
        if (not correct_path_result.is_initialized()) {
            INFO("No correct path");
            return false;
        }
//        set<scaffold_graph::ScaffoldVertex> permitted_vertices;
//        for (const auto &vertex: correct_path_result.get()) {
//            permitted_vertices.insert(vertex);
//            permitted_vertices.insert(vertex.getConjugateFromGraph(gp_.g));
//        }
//        permitted_vertices.erase(sink.getConjugateFromGraph(gp_.g));
//        permitted_vertices.erase(source.getConjugateFromGraph(gp_.g));
//        for (const auto &vertex: graph) {
//            if (permitted_vertices.find(vertex) == permitted_vertices.end()) {
//                INFO("Vertex " << vertex.int_id() << " is not permitted!");
//                return false;
//            }
//        }
    }

    return true;
}
PathClusterStatisticsExtractor::IdMap PathClusterStatisticsExtractor::GetIdMap(
        const PathClusterStatisticsExtractor::SimpleTransitionGraph &graph,
        const vector<PathClusterStatisticsExtractor::ScaffoldVertex> &correct_path) const {
    size_t current_id = 1;
    IdMap result;
//    if (not correct_path.empty()) {
//        for (const auto& vertex: correct_path) {
//            result.insert({vertex, std::to_string(current_id)});
//            result.insert({vertex.getConjugateFromGraph(gp_.g), std::to_string(current_id) + "'"});
//            ++current_id;
//        }
//        return result;
//    }
    for (const auto &vertex: graph) {
        result.insert({vertex, std::to_string(current_id)});
        ++current_id;
    }
    return result;
}
PathClusterStatisticsExtractor::ScaffoldEdgeMap PathClusterStatisticsExtractor::ConstructLengthMap(
        const PathClusterStatisticsExtractor::SimpleTransitionGraph &transition_graph,
        const PathClusterStatisticsExtractor::ScaffoldGraph &graph) const {
    PathClusterStatisticsExtractor::ScaffoldEdgeMap result;
    for (const auto &vertex: transition_graph) {
        if (result.find(vertex) == result.end()) {
            std::map<ScaffoldVertex, size_t> empty;
            result[vertex] = empty;
        }
        for (const auto &next: graph.OutgoingEdges(vertex)) {
            if (transition_graph.ContainsVertex(next.getEnd()) and
                transition_graph.ContainsEdge(vertex, next.getEnd())) {
                result.at(vertex)[next.getEnd()] = next.getLength();
            }
        }
    }
    return result;
}
SubgraphInfo::SubgraphInfo(const SubgraphInfo::SimpleTransitionGraph &graph,
                           const SubgraphInfo::ScaffoldVertex &source,
                           const SubgraphInfo::ScaffoldVertex &sink,
                           const PathClusterStorage &path_cluster_to_weight,
                           const vector<SubgraphInfo::VertexSet> &final_clusters,
                           const vector<vector<ScaffoldVertex>> &resulting_paths,
                           const vector<vector<ScaffoldVertex>> &all_paths,
                           const vector<SubgraphInfo::ScaffoldVertex> &correct_path,
                           const std::map<SubgraphInfo::ScaffoldVertex, string> &id_map,
                           const std::map<ScaffoldVertex, double> &vertex_to_cov,
                           const std::map<ScaffoldVertex, size_t> &vertex_to_len,
                           const ScaffoldEdgeMap &scaffold_edge_to_dist)
    : graph_(graph), source_(source), sink_(sink),
      path_cluster_to_weight_(path_cluster_to_weight),
      final_clusters_(final_clusters),
      resulting_paths_(resulting_paths),
      all_paths_(all_paths),
      correct_path_(correct_path), id_map_(id_map),
      vertex_to_cov_(vertex_to_cov),
      vertex_to_len_(vertex_to_len),
      scaffold_edge_to_dist_(scaffold_edge_to_dist) {}

ostream &operator<<(ostream &os, const SubgraphInfo &info) {
    const auto &graph = info.graph_;
    os << "Source: " << info.id_map_.at(info.source_) << "\n";
    os << "Sink: " << info.id_map_.at(info.sink_) << "\n";
    os << "Graph: " << "\n";
    for (const auto& vertex: graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            string current_short_id = info.id_map_.at(vertex);
            string next_short_id = info.id_map_.at(*it);
            os << current_short_id << " -> " << next_short_id << ", len: "
               << info.scaffold_edge_to_dist_.at(vertex).at(*it) << "\n";
        }
    }
    os <<"\nCorrect path\n";
    for (const auto& vertex: info.correct_path_) {
        os << info.id_map_.at(vertex) << " -> ";
    }
    os <<"\nVertex info:\n";
    for (const auto& vertex: info.graph_) {
        os << info.id_map_.at(vertex) << ", len: " << info.vertex_to_len_.at(vertex)
           << ", coverage: " << info.vertex_to_cov_.at(vertex) << "\n";
    }

    os <<"\nNormalized path clusters\n";
    for (const auto& entry: info.path_cluster_to_weight_) {
        const auto& vertex_set = entry.first;
        os << "{";
        for (const auto& v: vertex_set) {
            os << info.id_map_.at(v) << ", ";
        }
        os << "}: " << entry.second << "\n";
    }
    os << "\nFinal path clusters\n";
    for (const auto& cluster: info.final_clusters_) {
        os << "{";
        for (const auto& v: cluster) {
            os << info.id_map_.at(v) << ", ";
        }
        os << "}\n";
    }

    os <<"\nAll paths\n";
    for (const auto& path: info.all_paths_) {
        for (const auto &vertex: path) {
            os << info.id_map_.at(vertex) << " -> ";
        }
        os << "\n";
    }

    os <<"\nResulting paths\n";
    for (const auto& path: info.resulting_paths_) {
        for (const auto &vertex: path) {
            os << info.id_map_.at(vertex) << " -> ";
        }
        os << "\n";
    }
    return os;
}
void SubgraphInfoPrinter::PrintSubgraphInfo(const vector<SubgraphInfo> &info_collection, const string &output_path) const {
    std::ofstream fout(fs::append_path(output_path, "cluster_statistics"));
    for (const auto& info: info_collection) {
        fout << "----\n";
        fout << info;
        fout << "----\n";
    }
}
PathClusterExtractionParams::PathClusterExtractionParams(
        const Graph &g,
        shared_ptr<cluster_storage::InitialClusterStorage> init_cluster_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t linkage_distance)
    : g_(g),
      init_cluster_storage_(init_cluster_storage),
      barcode_extractor_(barcode_extractor),
      linkage_distance_(linkage_distance) {}
}