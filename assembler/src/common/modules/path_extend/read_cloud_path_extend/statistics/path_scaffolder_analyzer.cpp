//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_scaffolder_analyzer.hpp"

#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"

namespace path_extend {
namespace read_cloud {

boost::optional<size_t> PathDistanceEstimator::GetLongEdgeDistance(const PathDistanceEstimator::ScaffoldVertex &first,
                                                                   const PathDistanceEstimator::ScaffoldVertex &second) const {
    boost::optional<size_t> result;
    auto unique_predicate = [this](const EdgeId &edge) {
      return this->g_.length(edge) >= this->long_edge_threshold_ and this->reference_path_index_.Contains(edge);
    };
    auto last_unique_in_first = first.GetLastEdgeWithPredicate(unique_predicate);
    auto first_unique_in_second = second.GetFirstEdgeWithPredicate(unique_predicate);
    if (not first_unique_in_second.is_initialized() or not last_unique_in_first.is_initialized()) {
        DEBUG("No unique edges");
        return result;
    }
    auto first_index_entry = reference_path_index_.at(last_unique_in_first.get());
    auto second_index_entry = reference_path_index_.at(first_unique_in_second.get());
    size_t first_path = first_index_entry.path_;
    size_t second_path = second_index_entry.path_;
    if (first_path != second_path) {
        DEBUG("Different paths");
        return result;
    }
    size_t first_start = first_index_entry.start_pos_;
    size_t first_end = first_index_entry.end_pos_;
    size_t second_start = second_index_entry.start_pos_;
    size_t second_end = second_index_entry.end_pos_;
    if (second_start >= first_end) {
        result = second_start - first_end;
    } else {
        DEBUG("First edge: " << "[" << first_start << ", " << first_end << "]");
        DEBUG("Second edge: " << "[" << second_start << ", " << second_end << "]");
    }
    return result;
}
PathDistanceEstimator::PathDistanceEstimator(const Graph &g,
                                             const validation::ReferencePathIndex &reference_path_index,
                                             size_t long_edge_threshold)
    : g_(g), reference_path_index_(reference_path_index), long_edge_threshold_(long_edge_threshold) {}
PathPairInfo::PathPairInfo(size_t first_length, size_t second_length, size_t long_edge_distance, double score) :
    first_length_(first_length),
    second_length_(second_length),
    long_edge_distance_(long_edge_distance),
    score_(score) {}

std::ostream &operator<<(std::ostream &os, const PathPairInfo &info) {
    os << info.first_length_ << info.second_length_ << info.long_edge_distance_ << info.score_;
    return os;
}
PathPairDataset::PathPairDataset() {}
std::ostream &operator<<(std::ostream &os, const PathPairDataset &dataset) {
    os << "Length1" << "Length2" << "LongEdgeDist" << "Score" << std::endl;
    for (const auto &entry: dataset.data_) {
        os << entry << std::endl;
    }
    return os;
}
void PathPairDataset::Insert(const PathPairInfo &path_pair_info) {
    data_.push_back(path_pair_info);
}
PathPairDataset PathScaffolderAnalyzer::GetFalseNegativeDataset(const PathContainer &paths) const {
    const string path_to_reference = path_to_reference_;
    validation::ContigPathBuilder contig_path_builder(g_, index_, kmer_mapper_);

    auto named_reference_paths = contig_path_builder.GetContigPaths(path_to_reference);
    auto raw_reference_paths = contig_path_builder.StripNames(named_reference_paths);
    validation::FilteredReferencePathHelper path_helper(g_, index_, kmer_mapper_);
    auto filtered_reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference,
                                                                                    long_edge_length_threshold_);
    validation::GeneralTransitionStorageBuilder reference_transition_builder(g_, 1, false, false);
    auto reference_transitions = reference_transition_builder.GetTransitionStorage(filtered_reference_paths);
    INFO(reference_transitions.GetCoveredEdges().size() << " covered edges");
    INFO(reference_transitions.size() << " reference transitions");

    auto scaffold_vertices = ConstructScaffoldVertices(paths, reference_transitions);
    auto index = ConstructIndex(scaffold_vertices);
    auto graph = ConstructScaffoldGraph(scaffold_vertices, index);
    INFO("Graph contains " << graph.VertexCount() << " vertices and " << graph.EdgeCount() << " edges");
    const double score_threshold = 0.1;
    auto false_negative_edges = GetFalseNegativeEdges(graph, reference_transitions, score_threshold);
    INFO("False negative edges: " << false_negative_edges.size());

    validation::ReferencePathIndexBuilder path_index_builder;
    auto reference_path_index = path_index_builder.BuildReferencePathIndex(filtered_reference_paths);
    PathDistanceEstimator distance_estimator(g_, reference_path_index, long_edge_length_threshold_);
//    const size_t distance_threshold = 5000;

    PathPairDataset result;

    for (const auto &edge: false_negative_edges) {
        boost::optional<size_t>
            distance_result = distance_estimator.GetLongEdgeDistance(edge.getStart(), edge.getEnd());
        if (distance_result.is_initialized()) {
            INFO("Distance: " << distance_result.get());
            INFO("Score: " << edge.getWeight());
        } else {
            DEBUG("Distance undefined");
        }
    }
    return result;

}
std::shared_ptr<PathScaffolderAnalyzer::BarcodeIndex> PathScaffolderAnalyzer::ConstructIndex(
        const std::set<PathScaffolderAnalyzer::ScaffoldVertex> &vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper_, g_);
    auto graph_construction_params = configs_.scaff_con;
    const size_t tail_threshold = long_edge_length_threshold_;
    const size_t length_threshold = graph_construction_params.min_edge_length_for_barcode_collection;
    const size_t count_threshold = graph_construction_params.count_threshold;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    return scaffold_index_extractor;
}
PathScaffolderAnalyzer::ScaffoldGraph PathScaffolderAnalyzer::ConstructScaffoldGraph(
        const std::set<PathScaffolderAnalyzer::ScaffoldVertex> &vertices,
        std::shared_ptr<PathScaffolderAnalyzer::BarcodeIndex> index) const {
    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(g_, index);

    const double score_threshold = 0;
    INFO("Setting containment index threshold to " << score_threshold);

    auto initial_constructor = std::make_shared<scaffolder::ScoreFunctionScaffoldGraphConstructor>(g_,
                                                                                                   vertices,
                                                                                                   score_function,
                                                                                                   score_threshold,
                                                                                                   max_threads_);
    return *(initial_constructor->Construct());
}
std::set<PathScaffolderAnalyzer::ScaffoldVertex> PathScaffolderAnalyzer::ConstructScaffoldVertices(
        const PathContainer &paths,
        const validation::ContigTransitionStorage &transitions) const {
    std::set<ScaffoldVertex> path_set;
    for (const auto &path_pair: paths) {
        auto unique_predicate = [&transitions](const EdgeId &edge) {
          return transitions.IsEdgeCovered(edge);
        };
        if (std::any_of(path_pair.first->begin(), path_pair.first->end(), unique_predicate) or
            std::any_of(path_pair.second->begin(), path_pair.second->end(), unique_predicate)) {
            ScaffoldVertex first_vertex(path_pair.first);
            ScaffoldVertex second_vertex(path_pair.second);
            path_set.insert(first_vertex);
            path_set.insert(second_vertex);
            INFO("Length: " << first_vertex.GetLengthFromGraph(g_));
        }
    }
    return path_set;
}
std::set<PathScaffolderAnalyzer::ScaffoldEdge> PathScaffolderAnalyzer::GetFalseNegativeEdges(
        const PathScaffolderAnalyzer::ScaffoldGraph &graph,
        const validation::ContigTransitionStorage &transitions,
        double score_threshold) const {
    validation::ScaffoldGraphValidator validator(g_);
    auto is_covered = [&transitions](const EdgeId &edge) {
      return transitions.IsEdgeCovered(edge);
    };
    std::set<ScaffoldEdge> false_negative_edges;
    for (const ScaffoldEdge &edge: graph.edges()) {
        auto first_edge_result = edge.getStart().GetLastEdgeWithPredicate(is_covered);
        auto second_edge_result = edge.getEnd().GetFirstEdgeWithPredicate(is_covered);
        if (first_edge_result.is_initialized() and second_edge_result.is_initialized()) {
            transitions::Transition transition(first_edge_result.get(), second_edge_result.get());
            if (transitions.CheckTransition(transition) and math::le(edge.getWeight(), score_threshold)) {
                false_negative_edges.insert(edge);
            }
        }
    }
    return false_negative_edges;
}
PathScaffolderAnalyzer::PathScaffolderAnalyzer(const Graph &g,
                                               const debruijn_graph::Index &index,
                                               const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                               const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                               const ReadCloudConfigs &configs,
                                               const std::string &path_to_reference,
                                               size_t long_edge_length_threshold,
                                               size_t max_threads)
    : g_(g),
      index_(index),
      kmer_mapper_(kmer_mapper),
      barcode_mapper_(barcode_mapper),
      configs_(configs),
      path_to_reference_(path_to_reference),
      long_edge_length_threshold_(long_edge_length_threshold),
      max_threads_(max_threads) {}
}
}