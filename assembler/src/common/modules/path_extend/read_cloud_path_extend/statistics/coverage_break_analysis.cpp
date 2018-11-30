#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "coverage_break_analysis.hpp"
namespace path_extend {

GraphBreakAnalyzer::TransitionToBreaks GraphBreakAnalyzer::GetGraphBreaks(const string &path_to_reference,
                                                                          size_t length_threshold) const {
    auto raw_paths = GetRawPaths(path_to_reference);
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto long_edge_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);
    std::unordered_set<EdgeId> long_edges;
    for (const auto &path: long_edge_paths) {
        for (const auto &edge: path) {
            long_edges.insert(edge.edge_);
        }
    }
    validation::GeneralTransitionStorageBuilder reference_transition_builder(gp_.g, 1, false, false);
    auto reference_transitions = reference_transition_builder.BuildStorage(long_edge_paths);

    const double min_mapped_ratio = 0.5;
    auto fixed_long_edge_paths = FixLongEdgeAlignments(raw_paths, long_edges, min_mapped_ratio);

    validation::ReferencePathIndexBuilder path_index_builder;
    auto raw_path_index = path_index_builder.BuildReferencePathIndex(raw_paths);

    const size_t min_intersection = 1;
    size_t total_breaks = 0;
    std::unordered_map<transitions::Transition, std::vector<GraphBreak>> transition_to_breaks;
    for (const auto &transition: reference_transitions) {
        auto subpath = GetBarcodedSubpath(transition.first_, transition.second_, fixed_long_edge_paths,
                                          raw_path_index, min_intersection);
        auto graph_breaks = GetGraphBreaks(subpath);
        if (not graph_breaks.empty()) {
            transition_to_breaks.insert({transition, graph_breaks});
            total_breaks += graph_breaks.size();
        }
    }
    INFO(transition_to_breaks.size() << " broken transitions, " << total_breaks << " total breaks");
    return transition_to_breaks;
}
GraphBreakAnalyzer::GraphBreakAnalyzer(const conj_graph_pack &gp) : gp_(gp) {}
std::vector<GraphBreakAnalyzer::MappedPath> GraphBreakAnalyzer::FixLongEdgeAlignments(
        const std::vector<GraphBreakAnalyzer::MappedPath> &raw_paths,
        const std::unordered_set<EdgeId> &long_edges,
        double min_mapped_ratio) const {
    std::vector<MappedPath> result;
    for (const auto &path: raw_paths) {
        MappedPath fixed_path;
        for (const auto &ewm: path) {
            if (long_edges.find(ewm.edge_) != long_edges.end()) {
                auto mapped_length = static_cast<double>(ewm.mapping_.size());
                auto edge_length = static_cast<double>(gp_.g.length(ewm.edge_));
                VERIFY_DEV(not math::eq(edge_length, 0.0));
                double mapped_ratio = mapped_length / edge_length;
                if (math::ge(mapped_ratio, min_mapped_ratio)) {
                    fixed_path.push_back(ewm);
                }
            } else {
                fixed_path.push_back(ewm);
            }
        }
        result.push_back(fixed_path);
    }
    return result;
}
GraphBreakAnalyzer::SimplePath GraphBreakAnalyzer::GetBarcodedSubpath(
        const EdgeId &first, const EdgeId &second, const std::vector<GraphBreakAnalyzer::MappedPath> &raw_paths,
        const validation::ReferencePathIndex &ref_path_index, size_t min_intersection) const {
    SimplePath result;
    size_t path_id = ref_path_index.at(first).path_;
    VERIFY_DEV(ref_path_index.at(second).path_ == path_id);
    size_t first_pos = ref_path_index.at(first).edge_pos_;
    size_t second_pos = ref_path_index.at(second).edge_pos_;
    if (first_pos >= second_pos) {
        WARN("Backward reference transition");
        return result;
    }
    const auto &path = raw_paths[path_id];
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    result.push_back(first);
    for (size_t i = first_pos + 1; i < second_pos; ++i) {
        EdgeId middle_edge = path[i].edge_;
        size_t first_intersection = barcode_extractor->GetNumberOfSharedBarcodes(first, middle_edge);
        size_t second_intersection = barcode_extractor->GetNumberOfSharedBarcodes(middle_edge, second);
        if (std::min(first_intersection, second_intersection) >= min_intersection) {
            result.push_back(middle_edge);
        }
    }
    result.push_back(second);
    return result;
}
vector<GraphBreakAnalyzer::GraphBreak> GraphBreakAnalyzer::GetGraphBreaks(const GraphBreakAnalyzer::SimplePath &path) const {
    vector<GraphBreak> result;
    for (size_t i1 = 0, i2 = i1 + 1; i2 < path.size(); ++i1, ++i2) {
        VertexId first_end = gp_.g.EdgeEnd(path[i1]);
        VertexId second_start = gp_.g.EdgeStart(path[i2]);
        if (first_end != second_start) {
            result.emplace_back(path[i1], path[i2], first_end, second_start);
        }
    }
    return result;
}
void GraphBreakAnalyzer::PrintGraphBreaks(const GraphBreakAnalyzer::TransitionToBreaks &transition_to_breaks,
                                          const string &output_path) const {
    ofstream fout(output_path);
    for (const auto &entry: transition_to_breaks) {
        const auto &transition = entry.first;
        fout << transition.first_.int_id() << " " << transition.second_.int_id() << std::endl;
        for (const auto &graph_break: entry.second) {
            fout << graph_break.first_edge_ << " " << graph_break.second_edge_ << " " <<
                 graph_break.first_end_ << " " << graph_break.second_start_ << std::endl;
        }
    }
}
vector<GraphBreakAnalyzer::MappedPath> GraphBreakAnalyzer::GetRawPaths(const string &path_to_reference) const {
    validation::ContigPathBuilder contig_path_builder(gp_);
    auto named_raw_paths = contig_path_builder.GetRawPaths(path_to_reference);
    vector<vector<validation::EdgeWithMapping>> raw_paths;
    for (const auto &path: named_raw_paths) {
        vector<validation::EdgeWithMapping> raw_path;
        for (const auto &entry: path.mapping_path) {
            raw_path.emplace_back(entry.first, entry.second.mapped_range);
        }
        raw_paths.push_back(raw_path);
    }
    return raw_paths;
}
GraphBreakAnalyzer::GraphBreak::GraphBreak(const EdgeId &first_edge,
                                              const EdgeId &second_edge,
                                              const VertexId &first_end,
                                              const VertexId &second_start)
    : first_edge_(first_edge), second_edge_(second_edge), first_end_(first_end), second_start_(second_start) {}
}
