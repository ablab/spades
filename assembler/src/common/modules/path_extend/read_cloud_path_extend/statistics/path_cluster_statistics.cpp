#include "path_cluster_statistics.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"

namespace path_extend {

vector<SubgraphInfo> PathClusterStatisticsExtractor::GetAllSubgraphInfo(const ScaffoldGraphStorage &storage) {
    VERIFY_DEV(cfg::get().ts_res.debug_mode);
    vector<SubgraphInfo> result;
    ScaffoldGraphGapCloserParamsConstructor params_constructor;
    auto subgraph_extractor_params =
        params_constructor.ConstructSubgraphExtractorParamsFromConfig(storage.GetLargeLengthThreshold());
    auto path_extractor_params = params_constructor.ConstructPathClusterPredicateParamsFromConfig();
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
    auto cluster_storage_builder =
        std::make_shared<cluster_storage::EdgeInitialClusterStorageBuilder>(gp_.g, barcode_extractor_ptr,
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

    for (const auto& edge: univocal_edges) {
        auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(storage.GetSmallScaffoldGraph(),
                                                                          edge.getStart(), edge.getEnd());
        if (CheckSubgraph(subgraph, edge.getStart(), edge.getEnd(), validator)) {
            auto info =
                GetSubgraphInfo(subgraph, edge.getStart(), edge.getEnd(), validator, path_cluster_extraction_params);
            result.push_back(info);
        }
    }
    INFO("Analyzed " << result.size() << " subgraphs");
    return result;
}

SubgraphInfo PathClusterStatisticsExtractor::GetSubgraphInfo(
        const PathClusterStatisticsExtractor::SimpleTransitionGraph &graph,
        const scaffold_graph::ScaffoldVertex &source, const scaffold_graph::ScaffoldVertex &sink,
        const validation::SimpleTransitionGraphValidator &validator,
        const PathClusterStatisticsExtractor::PathClusterExtractionParams &path_cluster_extraction_params) const {
    PathClusterExtractorHelper path_cluster_extractor_helper(std::get<0>(path_cluster_extraction_params),
                                                             std::get<1>(path_cluster_extraction_params),
                                                             std::get<2>(path_cluster_extraction_params),
                                                             std::get<3>(path_cluster_extraction_params));
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
    DEBUG(path_clusters.size() << " path clusters");

    std::map<SubgraphInfo::VertexSet, size_t> cluster_to_weight;

    for (const auto& cluster: path_clusters) {
        const auto vertex_set = cluster.GetVertexSet();
        if (cluster_to_weight.find(vertex_set) == cluster_to_weight.end()) {
            cluster_to_weight[vertex_set] = 1;
        } else {
            cluster_to_weight[vertex_set]++;
        }
    }

    const auto correct_path = validator.GetCorrectPath(graph, source, sink);
    VERIFY_DEV(correct_path.is_initialized());

    SubgraphInfo result(graph, source, sink, cluster_to_weight, correct_path.get());
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
                                                   const validation::SimpleTransitionGraphValidator &validator) const {
    const size_t subgraph_size_threshold = 4;
    bool is_graph_large = graph.size() >= subgraph_size_threshold and graph.GetEdgesCount() >= graph.size();
    bool has_source_and_sink = graph.ContainsVertex(source) and graph.ContainsVertex(sink);
    if (not (is_graph_large and has_source_and_sink)) {
        return false;
    }
    bool has_correct_path = validator.GetCorrectPath(graph, source, sink).is_initialized();
    return is_graph_large and has_correct_path;
}
SubgraphInfo::SubgraphInfo(const SubgraphInfo::SimpleTransitionGraph &graph,
                           const SubgraphInfo::ScaffoldVertex &source,
                           const SubgraphInfo::ScaffoldVertex &sink,
                           const map<SubgraphInfo::VertexSet, size_t> &path_cluster_to_weight,
                           const vector<SubgraphInfo::ScaffoldVertex> &correct_path)
    : graph_(graph), source_(source), sink_(sink),
      path_cluster_to_weight_(path_cluster_to_weight), correct_path_(correct_path) {}

ostream &operator<<(ostream &os, const SubgraphInfo &info) {

    std::map<SubgraphInfo::ScaffoldVertex, size_t> vertex_to_short_id;
    const auto &graph = info.graph_;
    size_t current_id = 1;
    for (const auto &vertex: graph) {
        vertex_to_short_id.insert({vertex, current_id});
        ++current_id;
    }
    os << "Source: " << vertex_to_short_id.at(info.source_) << "\n";
    os << "Sink: " << vertex_to_short_id.at(info.sink_) << "\n";
    os << "Graph: " << "\n";
    for (const auto& vertex: graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            size_t current_short_id = vertex_to_short_id.at(vertex);
            size_t next_short_id = vertex_to_short_id.at(*it);
            os << current_short_id << " -> " << next_short_id << "\n";
        }
    }
    os <<"\nCorrect path\n";
    for (const auto& vertex: info.correct_path_) {
        os << vertex_to_short_id.at(vertex) << " -> ";
    }

    os <<"\nPath clusters\n";
    for (const auto& entry: info.path_cluster_to_weight_) {
        const auto& vertex_set = entry.first;
        os << "{";
        for (const auto& v: vertex_set) {
            os << vertex_to_short_id.at(v) << ", ";
        }
        os << "}: " << entry.second << "\n";
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
}