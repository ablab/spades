#include "reference_path_checker.hpp"

namespace cont_index {
void ReferencePathChecker::CheckAssemblyGraph(const std::string &reference_path) const {
    std::unordered_set<std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId>> correct_transitions;
    std::vector<std::vector<debruijn_graph::EdgeId>> reference_paths;
    std::ifstream ref_stream(reference_path);
    std::string path_string;
    std::string ref_name;
    size_t path_length;
    std::string edge_name;
    while (ref_stream >> path_string) {
        std::vector<debruijn_graph::EdgeId> reference_path;
        ref_name = path_string.substr(0, path_string.find(":"));
        INFO(ref_name);
        ref_stream >> path_length;
        INFO(path_length);
        for (size_t i = 0; i < path_length; ++i) {
            ref_stream >> edge_name;
//                INFO(edge_name);
            EdgeId edge = seg_to_edge_.at(edge_name);
//                INFO(edge.int_id());
            reference_path.push_back(edge);
        }
        INFO("Read reference path");
        reference_paths.push_back(reference_path);
    }
    INFO("Getting gap paths");
    size_t total_pairs = 0;
    size_t connected_pairs = 0;
    for (const auto &path: reference_paths) {
        if (path.size() < 2)
            continue;
        for (auto first = path.begin(), second = std::next(path.begin()); second != path.end(); ++first, ++second) {
            ++total_pairs;
            EdgeId first_edge = *first;
            EdgeId second_edge = *second;
            if (graph_.EdgeEnd(first_edge) == graph_.EdgeEnd(second_edge)) {
                ++connected_pairs;
            }
        }
    }
    INFO(connected_pairs << " connected pairs out of " << total_pairs);
}
}

