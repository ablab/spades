#include "common/pipeline/config_struct.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "split_index_statistics.hpp"

namespace path_extend {

SplitStatistics::SplitStatistics(const vector<SplitEntry> &data) : data_(data) {}
void SplitStatistics::Serialize(const string &path) {
    ofstream fout(path);
    fout << "SplitIndex,Status" << std::endl;
    for (const auto &entry: data_) {
        fout << entry.split_index_ << "," << entry.status_ << "\n";
    }
}
SplitStatistics SplitStatisticsExtractor::GetSplitStatistics(const string &path_to_reference,
                                                             size_t length_threshold) const {
    path_extend::validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);
    path_extend::validation::GeneralTransitionStorageBuilder forward_transition_builder(gp_.g, 1, false, false);
    auto reference_transition_storage = forward_transition_builder.GetTransitionStorage(reference_paths);
    std::unordered_set<Transition> reference_transitions;
    for (const auto &transition: reference_transition_storage) {
        reference_transitions.insert(transition);
    }

    path_extend::validation::GeneralTransitionStorageBuilder close_transition_builder(gp_.g, 5, false, false);
    auto close_transition_storage = close_transition_builder.GetTransitionStorage(reference_paths);
    std::unordered_set<Transition> close_transitions;
    for (const auto &transition: close_transition_storage) {
        if (reference_transitions.find(transition) == reference_transitions.end()) {
            close_transitions.insert(transition);
        }
    }

    path_extend::validation::GeneralTransitionStorageBuilder conj_transition_builder(gp_.g, 5, true, true);
    auto conj_transition_storage = conj_transition_builder.GetTransitionStorage(reference_paths);
    std::unordered_set<Transition> conj_transitions;
    for (const auto &transition: conj_transition_storage) {
        if (reference_transitions.find(transition) == reference_transitions.end() and
                close_transitions.find(transition) == close_transitions.end()) {
            conj_transitions.insert(transition);
        }
    }
    INFO(reference_transitions.size() << " reference transitions");
    INFO(close_transitions.size() << " close transitions");
    INFO(conj_transitions.size() << " conjugate transitions");

    std::set<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    for (const auto &path: reference_paths) {
        for (const auto &edge: path) {
            scaffold_vertices.insert(edge.edge_);
        }
    }
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const double EDGE_LENGTH_FRACTION = 0.5;
    const size_t count_threshold = 1;
    const size_t max_threads = cfg::get().max_threads;

    auto fraction_tail_threshold_getter =
        make_shared<barcode_index::FractionTailThresholdGetter>(gp_.g, EDGE_LENGTH_FRACTION);
    auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor,
                                                                           fraction_tail_threshold_getter,
                                                                           count_threshold, length_threshold,
                                                                           max_threads, scaffold_vertices);
    auto split_scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);

    vector<SplitEntry> data;

    for (const auto &transition: reference_transitions) {
        double split_index = GetSplitIndex(transition, split_scaffold_index_extractor);
        data.emplace_back(split_index, "correct");
    }

    for (const auto &transition: close_transitions) {
        double split_index = GetSplitIndex(transition, split_scaffold_index_extractor);
        data.emplace_back(split_index, "close");
    }

    for (const auto &transition: conj_transitions) {
        double split_index = GetSplitIndex(transition, split_scaffold_index_extractor);
        data.emplace_back(split_index, "conj");
    }

    SplitStatistics result(data);
    return result;
}
double SplitStatisticsExtractor::GetSplitIndex(const SplitStatisticsExtractor::Transition &transition,
                                               shared_ptr<SplitStatisticsExtractor::BarcodeExtractor> barcode_extractor) const {
    auto first = transition.first_;
    auto second = transition.second_;
    auto first_start = barcode_extractor->GetHeadEntry(first);
    auto first_end = barcode_extractor->GetTailEntry(first);
    auto second_start = barcode_extractor->GetHeadEntry(second);
    auto second_end = barcode_extractor->GetTailEntry(second);

    std::vector<barcode_index::BarcodeId> start_start_intersection;
    std::set_intersection(first_start.begin(), first_start.end(), second_start.begin(), second_start.end(),
                          std::back_inserter(start_start_intersection));
    std::vector<barcode_index::BarcodeId> start_end_intersection;
    std::set_intersection(first_start.begin(), first_start.end(), second_end.begin(), second_end.end(),
                          std::back_inserter(start_end_intersection));
    std::vector<barcode_index::BarcodeId> end_start_intersection;
    std::set_intersection(first_end.begin(), first_end.end(), second_start.begin(), second_start.end(),
                          std::back_inserter(end_start_intersection));
    std::vector<barcode_index::BarcodeId> end_end_intersection;
    std::set_intersection(first_end.begin(), first_end.end(), second_end.begin(), second_end.end(),
                          std::back_inserter(end_end_intersection));

    double split_index = 1.0;
    size_t max_false_intersection = std::max(start_start_intersection.size(),
                                             std::max(start_end_intersection.size(),
                                             end_end_intersection.size()));
    if (max_false_intersection != 0) {
        split_index = static_cast<double>(end_start_intersection.size()) / static_cast<double>(max_false_intersection);
    }
    return split_index;
}
void SplitStatisticsExtractor::ConstructAndSerialize(const string &path_to_reference,
                                                     const string &output_base,
                                                     size_t length_threshold) const {
    auto split_statistics = GetSplitStatistics(path_to_reference, length_threshold);
    const string output_path = fs::append_path(output_base, "split_statistics.csv");
    split_statistics.Serialize(output_path);
}
SplitStatisticsExtractor::SplitStatisticsExtractor(const conj_graph_pack &gp) : gp_(gp) {}
SplitEntry::SplitEntry(double split_index, const string &status) : split_index_(split_index), status_(status) {}
}