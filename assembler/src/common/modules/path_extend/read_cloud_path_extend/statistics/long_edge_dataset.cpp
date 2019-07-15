#include <random>
#include "long_edge_dataset.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
namespace path_extend {
namespace read_cloud {

void LongEdgePairDataset::Serialize(const string &path) {
    ofstream fout(path);
    fout <<
         "LeftId,LeftLength,LeftCov,LeftSize,RightId,RightLength,RightCov,RightSize,Intersection,Distance,Genome,Correct"
         << std::endl;
    for (const auto &entry: dataset_) {
        fout << entry.first_entry_.id_ << "," << entry.first_entry_.length_ << "," << entry.first_entry_.coverage_
             << "," << entry.first_entry_.barcodes_ << "," << entry.second_entry_.id_
             << "," << entry.second_entry_.length_ << "," << entry.second_entry_.coverage_
             << "," << entry.second_entry_.barcodes_ << "," << entry.intersection_ << "," << entry.distance_
             << "," << entry.genome_ << "," << entry.correct_ << std::endl;
    }
}
LongEdgePairDataset LongEdgePairDatasetExtractor::GetLongEdgeDataset(
    const vector<vector<validation::EdgeWithMapping>> &reference_paths) const {
    INFO("Getting long edge dataset")
    validation::ReferencePathIndexBuilder path_index_builder;
    auto reference_index = path_index_builder.BuildReferencePathIndex(reference_paths);
    validation::GeneralTransitionStorageBuilder forward_transition_builder(gp_.g, 1, false, false);
    auto reference_transition_storage = forward_transition_builder.GetTransitionStorage(reference_paths);
    const size_t close_distance = 5;
    validation::GeneralTransitionStorageBuilder close_transition_builder(gp_.g, close_distance, true, true);
    auto close_transition_storage = close_transition_builder.GetTransitionStorage(reference_paths);
    auto covered_edges_set = reference_transition_storage.GetCoveredEdges();

    auto distance_map = GetDistanceMap(reference_paths);
    vector<EdgeId> reference_edges(covered_edges_set.begin(), covered_edges_set.end());

    vector<LongEdgePairEntry> dataset;
    auto long_edge_extractor = ConstructLongEdgeExtractor();
    auto correct_entries = GetCorrectEntries(long_edge_extractor, reference_transition_storage,
                                             reference_index, distance_map);
    std::move(correct_entries.begin(), correct_entries.end(), std::back_inserter(dataset));

    INFO(reference_edges.size() << " reference edges");
    INFO(reference_paths.size() << " reference paths");
    INFO(dataset.size() << " correct pairs");

    const size_t sample_size = 100000;
    std::random_device rd;
    std::mt19937 generator(rd());
    size_t range_size = reference_edges.size();
    INFO("Covered edges: " << range_size);
    VERIFY(range_size > 0);
    std::uniform_int_distribution<size_t> distribution(0, range_size - 1);
    size_t counter = 0;
    size_t random_pairs = 0;
    while (counter <= sample_size) {
        size_t first_idx = distribution(generator);
        size_t second_idx = distribution(generator);
        EdgeId first_edge = reference_edges[first_idx];
        EdgeId second_edge = reference_edges[second_idx];
        if (AreNotClose(close_transition_storage, first_edge, second_edge)) {
            dataset.push_back(GetLongEdgePairEntry(long_edge_extractor, first_edge, second_edge, 1000000, 0, false));
        }
        ++counter;
        if (counter % (sample_size / 10) == 0) {
            INFO("Processed " << counter << " pairs out of " << sample_size);
            INFO(random_pairs << " random pairs.");
        }
    }
    LongEdgePairDataset result(dataset);
    return result;
}
map<LongEdgePairDatasetExtractor::Transition, size_t> LongEdgePairDatasetExtractor::GetDistanceMap(
    const vector<vector<LongEdgePairDatasetExtractor::EdgeWithMapping>> &reference_paths) const {
    std::map<transitions::Transition, size_t> result;
    for (const auto &path: reference_paths) {
        for (auto first = path.begin(), second = std::next(first); second != path.end(); ++first, ++second) {
            size_t first_end = (*first).mapping_.end_pos;
            size_t second_beginning = (*second).mapping_.start_pos;
            if (second_beginning >= first_end) {
                size_t distance = second_beginning - first_end;
                transitions::Transition t(first->edge_, second->edge_);
                result.insert({t, distance});
            }
        }
    }
    INFO(result.size() << " distances counted");
    return result;
}
vector<LongEdgePairEntry> LongEdgePairDatasetExtractor::GetCorrectEntries(
    shared_ptr<LongEdgePairDatasetExtractor::BarcodeExtractor> long_edge_extractor,
    const validation::ContigTransitionStorage &reference_transition_storage,
    const validation::ReferencePathIndex &long_edge_path_index,
    const map<LongEdgePairDatasetExtractor::Transition, size_t> &distance_map) const {
    vector<LongEdgePairEntry> correct_entries;

    for (const auto &transition: reference_transition_storage) {
        DEBUG("Getting path ids");
        VERIFY(reference_transition_storage.IsEdgeCovered(transition.first_));
        VERIFY(reference_transition_storage.IsEdgeCovered(transition.second_));
        VERIFY((long_edge_path_index.Contains(transition.first_)));
        VERIFY((long_edge_path_index.Contains(transition.second_)));
        size_t first_path_id = long_edge_path_index.at(transition.first_).path_;
        size_t second_path_id = long_edge_path_index.at(transition.second_).path_;
        DEBUG(first_path_id << ", " << second_path_id);
        if (first_path_id != second_path_id) {
            WARN("Correct transition from different paths!");
        } else {
            if (distance_map.find(transition) != distance_map.end()) {
                size_t distance = distance_map.at(transition);
                correct_entries.push_back(GetLongEdgePairEntry(long_edge_extractor,
                                                               transition.first_, transition.second_, distance,
                                                               first_path_id, true));
            }
        }
    }
    return correct_entries;
}
LongEdgePairEntry LongEdgePairDatasetExtractor::GetLongEdgePairEntry(
    shared_ptr<LongEdgePairDatasetExtractor::BarcodeExtractor> long_edge_extractor,
    const EdgeId &first, const EdgeId &second,
    size_t distance, size_t path_id, bool correct) const {
    const size_t tail_threshold = scaffold_graph_storage_.GetSmallLengthThreshold();
    barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor(gp_.barcode_mapper_ptr, gp_.g);

    LongEdgeEntry first_entry(first.int_id(), gp_.g.length(first), gp_.g.coverage(first),
                              long_edge_extractor->GetTailSize(first));
    LongEdgeEntry second_entry(second.int_id(), gp_.g.length(second), gp_.g.coverage(second),
                               long_edge_extractor->GetHeadSize(second));
    size_t intersection = long_edge_extractor->GetIntersectionSize(first, second);
    LongEdgePairEntry result(first_entry, second_entry, intersection, distance, path_id, correct);

    auto first_conj = gp_.g.conjugate(first);
    auto head_barcodes = barcode_extractor.GetBarcodesFromHead(first_conj, 1, tail_threshold);
    auto tail_barcodes = barcode_extractor.GetBarcodesFromHead(second, 1, tail_threshold);
    vector<barcode_index::BarcodeId> intersection_b;
    std::set_intersection(head_barcodes.begin(), head_barcodes.end(), tail_barcodes.begin(), tail_barcodes.end(),
                          std::back_inserter(intersection_b));
    VERIFY_DEV(first_entry.barcodes_ == head_barcodes.size());
    VERIFY_DEV(second_entry.barcodes_ == tail_barcodes.size());
    VERIFY_DEV(result.intersection_ == intersection_b.size());

    return result;
}
bool LongEdgePairDatasetExtractor::AreNotClose(const validation::ContigTransitionStorage &close_transition_storage,
                                               const EdgeId &first,
                                               const EdgeId &second) const {
    bool are_close = first == second or gp_.g.conjugate(first) == second or
        close_transition_storage.CheckTransition(first, second);
    return not are_close;
}
shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> LongEdgePairDatasetExtractor::ConstructLongEdgeExtractor() const {
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
LongEdgePairDatasetExtractor::LongEdgePairDatasetExtractor(const conj_graph_pack &gp,
                                                           const ScaffoldGraphStorage &scaffold_graph_storage,
                                                           size_t max_threads) :
    gp_(gp),
    scaffold_graph_storage_(scaffold_graph_storage),
    max_threads_(max_threads) {}
LongEdgePairDataset LongEdgePairDatasetExtractor::GetLongEdgeDataset(const scaffold_graph::ScaffoldGraph &graph,
                                                                     const string &path_to_reference) const {
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromGraph(path_to_reference, graph);
    return GetLongEdgeDataset(reference_paths);
}
void LongEdgePairDatasetExtractor::ConstructAndSerialize(const string &path_to_reference,
                                                         const string &output_base) const {
    const string reference_path = path_to_reference;
    size_t long_threshold = scaffold_graph_storage_.GetSmallLengthThreshold();
    size_t ultralong_threshold = scaffold_graph_storage_.GetLargeLengthThreshold();
    INFO(scaffold_graph_storage_.GetSmallScaffoldGraph().VertexCount() << " long edges");
    auto long_edge_dataset = GetLongEdgeDataset(scaffold_graph_storage_.GetSmallScaffoldGraph(), reference_path);
    INFO(scaffold_graph_storage_.GetLargeScaffoldGraph().VertexCount() << " ultralong edges");
    auto ultralong_edge_dataset = GetLongEdgeDataset(scaffold_graph_storage_.GetLargeScaffoldGraph(), reference_path);
    const string output_name = "long_edge_dataset_";
    const string long_output_path = fs::append_path(output_base, output_name + std::to_string(long_threshold));
    const string ultralong_output_path = fs::append_path(output_base, output_name +
        std::to_string(ultralong_threshold));
    long_edge_dataset.Serialize(long_output_path);
    ultralong_edge_dataset.Serialize(ultralong_output_path);
}
}
}