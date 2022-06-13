//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "split_index_statistics.hpp"

#include "barcode_index/scaffold_vertex_index_builder.hpp"

namespace path_extend {
namespace read_cloud {

SplitStatistics::SplitStatistics(const std::vector<SplitEntry> &data) : data_(data) {}
void SplitStatistics::Serialize(const std::string &path) {
    std::ofstream fout(path);
    fout << "SplitIndex,Status" << std::endl;
    for (const auto &entry: data_) {
        fout << entry.split_index_ << "," << entry.status_ << "\n";
    }
}
SplitStatistics SplitStatisticsExtractor::GetSplitStatistics(const std::filesystem::path &path_to_reference,
                                                             size_t length_threshold) const {
    validation::FilteredReferencePathHelper path_helper(g_, index_, kmer_mapper_);
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, length_threshold, 5000.0);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromUnique(path_to_reference, unique_storage);

    validation::GeneralTransitionStorageBuilder forward_transition_builder(g_, 1, false, false);
    auto reference_transition_storage = forward_transition_builder.GetTransitionStorage(reference_paths);
    std::unordered_set<Transition> reference_transitions;
    for (const auto &transition: reference_transition_storage) {
        reference_transitions.insert(transition);
    }
    validation::GeneralTransitionStorageBuilder close_transition_builder(g_, 5, false, false);
    auto close_transition_storage = close_transition_builder.GetTransitionStorage(reference_paths);
    std::unordered_set<Transition> close_transitions;
    for (const auto &transition: close_transition_storage) {
        if (reference_transitions.find(transition) == reference_transitions.end()) {
            close_transitions.insert(transition);
        }
    }
    validation::GeneralTransitionStorageBuilder conj_transition_builder(g_, 5, true, true);
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
    for (const auto &path: reference_paths.paths_) {
        for (const auto &edge: path) {
            scaffold_vertices.insert(edge.edge_);
        }
    }
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper_, g_);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const double EDGE_LENGTH_FRACTION = 0.5;
    const size_t count_threshold = 1;

    auto fraction_tail_threshold_getter =
        std::make_shared<barcode_index::FractionTailThresholdGetter>(g_, EDGE_LENGTH_FRACTION);
    auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor,
                                                                           fraction_tail_threshold_getter,
                                                                           count_threshold, length_threshold,
                                                                           max_threads_, scaffold_vertices);
    auto split_scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);

    std::vector<SplitEntry> data;
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
double SplitStatisticsExtractor::GetSplitIndex(
        const SplitStatisticsExtractor::Transition &transition,
        std::shared_ptr<SplitStatisticsExtractor::BarcodeExtractor> barcode_extractor) const {
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
void SplitStatisticsExtractor::ConstructAndSerialize(const std::filesystem::path &path_to_reference,
                                                     const std::filesystem::path &output_base,
                                                     size_t length_threshold) const {
    auto split_statistics = GetSplitStatistics(path_to_reference, length_threshold);
    const std::string output_path = output_base / "split_statistics.csv";
    split_statistics.Serialize(output_path);
}
SplitStatisticsExtractor::SplitStatisticsExtractor(const graph_pack::GraphPack &gp,
                                                   size_t max_threads) :
    gp_(gp),
    g_(gp_.get<Graph>()),
    index_(gp_.get<EdgeIndex<Graph>>()),
    kmer_mapper_(gp_.get<KmerMapper<Graph>>()),
    barcode_mapper_(gp_.get<barcode_index::FrameBarcodeIndex<Graph>>()),
    max_threads_(max_threads) {}
SplitEntry::SplitEntry(double split_index, const std::string &status) : split_index_(split_index), status_(status) {}
}
}