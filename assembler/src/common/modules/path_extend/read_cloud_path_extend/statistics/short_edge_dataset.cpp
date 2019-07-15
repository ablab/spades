#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "read_cloud_path_extend/validation/reference_path_index.hpp"
#include "short_edge_dataset.hpp"

namespace path_extend {
namespace read_cloud {

ShortEdgeDataset ShortEdgeDatasetExtractor::GetShortEdgeDataset(
    const vector<vector<ShortEdgeDatasetExtractor::EdgeWithMapping>> &reference_paths,
    const vector<vector<ShortEdgeDatasetExtractor::EdgeWithMapping>> &filtered_reference_paths) const {
    vector<vector<validation::EdgeWithMapping>> current_paths = reference_paths;
    INFO("Getting short_edge_dataset");
    unordered_set<EdgeId> long_edges;
    for (const auto &path: filtered_reference_paths) {
        for (const auto &ewm: path) {
            long_edges.insert(ewm.edge_);
        }
    }
    validation::ReferencePathIndexBuilder path_index_builder;
    auto reference_index = path_index_builder.BuildReferencePathIndexForSet(reference_paths, long_edges);
    validation::StrictTransitionStorageBuilder transition_storage_builder;
    auto transition_storage = transition_storage_builder.GetTransitionStorage(filtered_reference_paths);
    auto entries = GetShortEdgeEntries(transition_storage, reference_index, reference_paths);
    ShortEdgeDataset short_edge_dataset(entries);
    return short_edge_dataset;
}
vector<ShortEdgeEntry> ShortEdgeDatasetExtractor::GetShortEdgeEntries(
    const validation::ContigTransitionStorage &transition_storage,
    const validation::ReferencePathIndex &long_edge_path_index,
    const vector<vector<ShortEdgeDatasetExtractor::EdgeWithMapping>> &reference_paths) const {
    vector<ShortEdgeEntry> entries;
    INFO(transition_storage.size() << " correct long edge transitions");
    const size_t MAX_RANDOM_EDGES = 5000;
    size_t total_correct_edges = 0;
    size_t total_random_edges = 0;
    for (const auto &transition: transition_storage) {
        EdgeId first = transition.first_;
        EdgeId second = transition.second_;
        VERIFY(long_edge_path_index.Contains(first));
        VERIFY(long_edge_path_index.Contains(second));
        size_t first_path_id = long_edge_path_index.at(first).path_;
        size_t second_path_id = long_edge_path_index.at(second).path_;
        if (first_path_id != second_path_id) {
//                WARN("Reference transition edges" << first.int_id() << " and " << second.int_id()
//                                                  << "belong to different references. Skipping.");
            continue;
        }
        const vector<EdgeWithMapping> &reference_path = reference_paths[first_path_id];
        size_t first_pos = long_edge_path_index.at(first).edge_pos_;
        size_t second_pos = long_edge_path_index.at(second).edge_pos_;
        auto long_edge_barcode_extractor = ConstructLongEdgeExtractor();
        auto first_entry = long_edge_barcode_extractor->GetTailEntry(first);
        auto second_entry = long_edge_barcode_extractor->GetHeadEntry(second);
        auto correct_edges = GetEdgesBetweenPair(first_pos, second_pos, reference_path);
        auto random_edges = GetReachableEdges(first);
        total_correct_edges += correct_edges.size();
        total_random_edges += random_edges.size();
        DEBUG(correct_edges.size() << " correct edges.");
        DEBUG(random_edges.size() << " random edges.");
        for (const auto &edge: correct_edges) {
            auto short_edge_entry = GetShortEdgeEntry(edge, first_entry, second_entry,
                                                      gp_.g.coverage(first), gp_.g.coverage(second), true);
            entries.push_back(short_edge_entry);
        }
        size_t current_random_edges = 0;
        for (const auto &edge: random_edges) {
            if (correct_edges.find(edge) == correct_edges.end()) {
                auto short_edge_entry = GetShortEdgeEntry(edge, first_entry, second_entry,
                                                          gp_.g.coverage(first), gp_.g.coverage(second), false);
                entries.push_back(short_edge_entry);
                ++current_random_edges;
            }
            if (current_random_edges > MAX_RANDOM_EDGES) {
                break;
            }
        }
    }
    INFO("Total correct edges: " << total_correct_edges);
    INFO("Total random edges: " << total_random_edges);
    return entries;
}

ShortEdgeEntry ShortEdgeDatasetExtractor::GetShortEdgeEntry(EdgeId short_edge,
                                                            const barcode_index::SimpleVertexEntry &left_entry,
                                                            const barcode_index::SimpleVertexEntry &right_entry,
                                                            double left_coverage,
                                                            double right_coverage,
                                                            bool correct) const {
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr,
                                                                                        gp_.g);
    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(gp_.g, barcode_extractor);
    size_t length = gp_.g.length(short_edge);
    double coverage = gp_.g.coverage(short_edge);
    size_t left_intersection = short_edge_extractor->GetIntersectionSize(short_edge, left_entry);
    size_t right_intersection = short_edge_extractor->GetIntersectionSize(short_edge, right_entry);
    size_t barcodes = short_edge_extractor->GetHeadSize(short_edge);
    size_t left_size = left_entry.size();
    size_t right_size = right_entry.size();
    ShortEdgeEntry entry(short_edge.int_id(), left_size, right_size, barcodes, left_intersection,
                         right_intersection, left_coverage, right_coverage, length, coverage, correct);
    return entry;
}

std::unordered_set<EdgeId> ShortEdgeDatasetExtractor::GetReachableEdges(const EdgeId &long_edge) const {
    const size_t DISTANCE_BOUND = 40000;
    size_t min_length = scaffold_graph_storage_.GetSmallLengthThreshold();
    unordered_set<EdgeId> reached_edges;
    DijkstraHelper<Graph> helper;
    auto unique_dijkstra = helper.CreateLengthBoundedDijkstra(gp_.g, DISTANCE_BOUND, min_length);
    unique_dijkstra.Run(gp_.g.EdgeEnd(long_edge));
    for (const auto &reached_vertex: unique_dijkstra.ReachedVertices()) {
        const auto &outgoing_edges = gp_.g.OutgoingEdges(reached_vertex);
        for (const auto &edge: outgoing_edges) {
            reached_edges.insert(edge);
        }
    }
    return reached_edges;
}

std::unordered_set<EdgeId> ShortEdgeDatasetExtractor::GetEdgesBetweenPair(size_t first_pos,
                                                                          size_t second_pos,
                                                                          const vector<EdgeWithMapping> &reference_path) const {
    std::unordered_set<EdgeId> correct_edges;
    for (size_t i = first_pos + 1; i < second_pos; ++i) {
        EdgeId middle = reference_path[i].edge_;
        correct_edges.insert(middle);
    }
    return correct_edges;
}
shared_ptr<ShortEdgeDatasetExtractor::BarcodeExtractor> ShortEdgeDatasetExtractor::ConstructLongEdgeExtractor() const {
    size_t min_length = scaffold_graph_storage_.GetSmallLengthThreshold();
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr,
                                                                                        gp_.g);
    const size_t tail_threshold = min_length;
    const size_t length_threshold = 500;
    const size_t count_threshold = 1;
    auto tail_threshold_getter = make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    const auto &scaffold_graph = scaffold_graph_storage_.GetSmallScaffoldGraph();
    set<scaffold_graph::ScaffoldVertex> vertices;
    for (const auto &vertex: scaffold_graph.vertices()) {
        vertices.insert(vertex);
    }
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, vertices);
    auto scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    return scaffold_index_extractor;
}

ShortEdgeDataset ShortEdgeDatasetExtractor::GetShortEdgeDataset(size_t length_threshold,
                                                                const string &path_to_reference) const {
    validation::ContigPathBuilder contig_path_builder(gp_);
    auto named_reference_paths = contig_path_builder.GetContigPaths(path_to_reference);
    auto reference_paths = contig_path_builder.StripNames(named_reference_paths);
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto
        filtered_reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);
    return GetShortEdgeDataset(reference_paths, filtered_reference_paths);
}
void ShortEdgeDatasetExtractor::ConstructAndSerialize(const string &path_to_reference,
                                                      const string &output_base) const {
    const string reference_path = path_to_reference;
    size_t long_threshold = scaffold_graph_storage_.GetSmallLengthThreshold();
    size_t ultralong_threshold = scaffold_graph_storage_.GetLargeLengthThreshold();
    INFO(scaffold_graph_storage_.GetSmallScaffoldGraph().VertexCount() << " long edges");
    auto short_long_edge_dataset = GetShortEdgeDataset(long_threshold, reference_path);
    INFO(scaffold_graph_storage_.GetLargeScaffoldGraph().VertexCount() << " ultralong edges");
    auto short_ultralong_edge_dataset = GetShortEdgeDataset(ultralong_threshold, reference_path);
    const string output_name = "short_edge_dataset_";
    const string long_output_path = fs::append_path(output_base, output_name + std::to_string(long_threshold));
    const string ultralong_output_path = fs::append_path(output_base, output_name +
        std::to_string(ultralong_threshold));
    short_long_edge_dataset.Serialize(long_output_path);
    short_ultralong_edge_dataset.Serialize(ultralong_output_path);
}
}
}