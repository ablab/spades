//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "transition_extractor.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

std::vector<NamedSimplePath> ContigPathBuilder::GetContigPaths(const string &path_to_contigs) const {
    auto raw_contig_paths = GetRawPaths(path_to_contigs);
    DEBUG(raw_contig_paths.size() << " raw paths");
    std::vector<NamedSimplePath> fixed_paths = FixMappingPaths(raw_contig_paths);
    DEBUG(fixed_paths.size() << " fixed paths");
    return fixed_paths;
}
std::vector<std::vector<EdgeWithMapping>> ContigPathBuilder::StripNames(
        const std::vector<NamedSimplePath> &named_paths) const {
    std::vector<std::vector<EdgeWithMapping>> result;
    for (const auto &named_path: named_paths) {
        result.push_back(named_path.path_);
    }
    return result;
}
std::vector<NamedPath> ContigPathBuilder::GetRawPaths(const string &contig_path) const {
    std::vector<NamedPath> contig_paths;
    const debruijn_graph::Index &index = gp_.index;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper = gp_.kmer_mapper;
    auto contig_stream_ptr = std::make_shared<io::FileReadStream>(contig_path);
    auto rc_contig_stream_ptr = io::RCWrap<io::SingleRead>(contig_stream_ptr);

    std::vector<io::SingleRead> contigs;
    const size_t contig_length_threshold = 4000;
    while (!rc_contig_stream_ptr->eof()) {
        io::SingleRead contig;
        (*rc_contig_stream_ptr) >> contig;
        auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::Index>>
            (gp_.g, index, kmer_mapper);
        MappingPath<EdgeId> path = mapper->MapRead(contig);
        if (path.size() > 0 and path.back().second.initial_range.end_pos > contig_length_threshold) {
            NamedPath named_path(path, RemoveSpacesFromName(contig.name()));
            contig_paths.push_back(named_path);
        }
    }
    return contig_paths;
}
string ContigPathBuilder::RemoveSpacesFromName(const string &name) const {
    string result;
    const char old_space = ' ';
    const char new_space = '_';
    for (const char c: name) {
        if (c != old_space) {
            result.push_back(c);
        } else {
            result.push_back(new_space);
        }
    }
    return result;
}
std::vector<NamedSimplePath> ContigPathBuilder::FixMappingPaths(const std::vector<NamedPath> &contig_paths) const {
    std::vector<NamedSimplePath> result;
    debruijn_graph::GappedPathExtractor gapped_path_extractor(gp_.g);
    for (const auto &named_path: contig_paths) {
        auto path = named_path.mapping_path;
        std::vector<std::vector<EdgeId>> fixed_paths = gapped_path_extractor(path);
        std::vector<EdgeWithMapping> long_path;
        size_t prefix_len = 0;
        for (const auto &fixed_path: fixed_paths) {
            for (const auto &edge: fixed_path) {
                Range mapping(prefix_len, prefix_len + gp_.g.length(edge));
                EdgeWithMapping ewm(edge, mapping);
                long_path.push_back(ewm);
                prefix_len += gp_.g.length(edge);
            }
        }
        NamedSimplePath named_simple_path(long_path, named_path.name);
        result.push_back(named_simple_path);
    }
    return result;
}
std::vector<std::vector<EdgeWithMapping>> ContigPathFilter::FilterPaths(const ReferencePaths &paths) const {
    ReferencePaths filtered_paths;
    for (const auto &path: paths) {
        std::vector<EdgeWithMapping> filtered_path;
        for (const auto &ewm: path) {
            if (edges_.find(ewm.edge_) != edges_.end()) {
                filtered_path.push_back(ewm);
            }
        }
        if (not filtered_path.empty()) {
            filtered_paths.push_back(filtered_path);
        }
    }
    auto merged_paths = MergeSameEdges(filtered_paths);

    auto result = RemoveRepeats(merged_paths);

    return result;
}
std::vector<std::vector<EdgeWithMapping>> ContigPathFilter::MergeSameEdges(const ReferencePaths &paths) const {
    ReferencePaths result;
    for (const auto &path: paths) {
        std::vector<EdgeWithMapping> current_path;
        for (const auto &ewm: path) {
            if (current_path.empty() or current_path.back().edge_ != ewm.edge_) {
                current_path.push_back(ewm);
            } else {
                if (not current_path.empty()) {
                    Range prev_mapping = current_path.back().mapping_;
                    Range new_mapping(prev_mapping.start_pos, ewm.mapping_.end_pos);
                    EdgeWithMapping new_ewm(ewm.edge_, new_mapping);
                    current_path.back() = new_ewm;
                }
            }
        }

        //merge last and first edge for cycle paths
        auto last_ewm = current_path.back();
        auto first_ewm = current_path[0];
        if (last_ewm.edge_ == first_ewm.edge_ and current_path.size() > 1) {
            Range new_mapping(std::min(last_ewm.mapping_.start_pos, first_ewm.mapping_.start_pos),
                              std::max(last_ewm.mapping_.end_pos, first_ewm.mapping_.end_pos));
            current_path[0].mapping_ = new_mapping;
            current_path.pop_back();
        }
        if (not current_path.empty()) {
            result.push_back(current_path);
        }
    }
    return result;
}
std::unordered_set<EdgeId> ContigPathFilter::GetRepeats(const ReferencePaths &paths) const {
    std::unordered_set<EdgeId> repeats;
    std::unordered_set<EdgeId> visited;
    for (const auto &path: paths) {
        for (const auto &edge: path) {
            bool is_visited = not(visited.insert(edge.edge_).second);
            if (is_visited) {
                repeats.insert(edge.edge_);
            }
        }
    }
    return repeats;
}
std::vector<std::vector<EdgeWithMapping>> ContigPathFilter::RemoveRepeats(const ReferencePaths &paths) const {
    std::vector<std::vector<EdgeWithMapping>> result;
    auto repeats = GetRepeats(paths);
    INFO("Detected " << repeats.size() << " repeats");
    for (const auto &path: paths) {
        std::vector<EdgeWithMapping> new_path;
        for (const auto &edge: path) {
            if (repeats.find(edge.edge_) == repeats.end()) {
                new_path.push_back(edge);
            }
        }
        result.push_back(new_path);
    }
    return result;
}

ContigTransitionStorage StrictTransitionStorageBuilder::BuildStorage(
        const std::vector<std::vector<EdgeWithMapping>> &long_edges) const {
    ContigTransitionStorage storage;
    for (const auto &path: long_edges) {
        for (auto it1 = path.begin(), it2 = std::next(path.begin());
             it1 != path.end() and it2 != path.end(); ++it1, ++it2) {
            EdgeId edge1 = (*it1).edge_;
            EdgeId edge2 = (*it2).edge_;
            storage.InsertEdge(edge1);
            storage.InsertEdge(edge2);
            storage.InsertTransition(edge1, edge2);
        }
    }
    return storage;
}
ContigTransitionStorage ReverseTransitionStorageBuilder::BuildStorage(
        const std::vector<std::vector<EdgeWithMapping>> &long_edges) const {
    ContigTransitionStorage storage;
    for (const auto &path: long_edges) {
        for (auto it1 = path.rbegin(), it2 = std::next(path.rbegin());
             it1 != path.rend() and it2 != path.rend(); ++it1, ++it2) {
            EdgeId edge1 = (*it1).edge_;
            EdgeId edge2 = (*it2).edge_;
            storage.InsertEdge(edge1);
            storage.InsertEdge(edge2);
            storage.InsertTransition(edge1, edge2);
        }
    }
    return storage;
}
ContigTransitionStorage ConjugateTransitionStorageBuilder::BuildStorage(
        const std::vector<std::vector<EdgeWithMapping>> &long_edges) const {
    ContigTransitionStorage storage;
    for (const auto &path: long_edges) {
        for (auto it1 = path.rbegin(), it2 = std::next(path.rbegin());
             it1 != path.rend() and it2 != path.rend(); ++it1, ++it2) {
            EdgeId edge1 = (*it1).edge_;
            EdgeId edge2 = g_.conjugate((*it2).edge_);
            storage.InsertEdge(edge1);
            storage.InsertEdge(edge2);
            storage.InsertTransition(edge1, edge2);
        }
        for (auto it1 = path.begin(), it2 = std::next(path.begin());
             it1 != path.end() and it2 != path.end(); ++it1, ++it2) {
            EdgeId edge1 = (*it1).edge_;
            EdgeId edge2 = g_.conjugate((*it2).edge_);
            storage.InsertEdge(edge1);
            storage.InsertEdge(edge2);
            storage.InsertTransition(edge1, edge2);
        }
    }
    return storage;
}
ContigTransitionStorage GeneralTransitionStorageBuilder::BuildStorage(
        const std::vector<std::vector<EdgeWithMapping>> &paths) const {
    ContigTransitionStorage storage;
    for (const auto &path: paths) {
        for (size_t prev_idx = 0; prev_idx < path.size(); ++prev_idx) {
            EdgeId first = path[prev_idx].edge_;
            storage.InsertEdge(first);
            for (size_t next_idx = prev_idx + 1; next_idx != path.size() and next_idx <= prev_idx + distance_;
                 ++next_idx) {
                EdgeId second = path[next_idx].edge_;
                ProcessPair(first, second, with_conjugate_, with_reverse_, storage);
            }
            if (prev_idx + distance_ >= path.size() and prev_idx + distance_ < 2 * path.size()) {
                size_t rem = prev_idx + distance_ - path.size();
                for (size_t right_idx = 0; right_idx <= rem; ++right_idx) {
                    EdgeId second = path[right_idx].edge_;
                    ProcessPair(first, second, with_conjugate_, with_reverse_, storage);
                }
            }
        }
    }
    return storage;
}
void GeneralTransitionStorageBuilder::ProcessPair(EdgeId first,
                                                  EdgeId second,
                                                  bool with_conjugate,
                                                  bool with_reverse,
                                                  ContigTransitionStorage &storage) const {
    storage.InsertEdge(second);
    storage.InsertTransition(first, second);
    if (with_conjugate) {
        EdgeId conjugate = g_.conjugate(second);
        storage.InsertEdge(conjugate);
        storage.InsertTransition(first, conjugate);
    }
    if (with_reverse) {
        storage.InsertTransition(second, first);
        if (with_conjugate) {
            storage.InsertEdge(g_.conjugate(first));
            storage.InsertTransition(second, g_.conjugate(first));
        }
    }

}
ContigTransitionStorage ApproximateTransitionStorageBuilder::BuildStorage(
        const std::vector<std::vector<EdgeWithMapping>> &long_edges) const {
    ContigTransitionStorage storage;
    const size_t threshold = 5000;
    const size_t mapping_error_threshold = 500;
    for (const auto &path: long_edges) {
        for (const auto &ewm: path) {
            storage.InsertEdge(ewm.edge_);
        }
    }
    for (const auto &path: long_edges) {
        for (auto it = path.begin(); it != path.end(); ++it) {
            EdgeId first = (*it).edge_;
            size_t pos = (*it).mapping_.end_pos;
            for (auto it_next = std::next(it); it_next != path.end(); ++it_next) {
                EdgeId second = (*it_next).edge_;
                size_t next_pos = (*it_next).mapping_.start_pos;
                VERIFY(next_pos + mapping_error_threshold >= pos);
                if (pos + threshold > next_pos) {
                    storage.InsertTransition(first, second);
                } else {
                    break;
                }
            }
        }
    }
    return storage;
}
std::vector<ClusterTransitionExtractor::Transition> ClusterTransitionExtractor::ExtractAllTransitionsFromNonPathCluster(
        const cluster_storage::Cluster &cluster) {
    const auto &internal_graph = cluster.GetInternalGraph();
    std::vector<Transition> result;
    for (const auto &vertex: internal_graph) {
        for (auto it = internal_graph.outcoming_begin(vertex); it != internal_graph.outcoming_end(vertex); ++it) {
            result.emplace_back(vertex, *it);
        }
    }
    return result;
}
std::vector<ClusterTransitionExtractor::Transition> ClusterTransitionExtractor::ExtractGoodTransitionsFromNonPathCluster(
        const cluster_storage::Cluster &cluster) {
    const auto &internal_graph = cluster.GetInternalGraph();
    std::vector<Transition> result;
    for (const auto &start: internal_graph) {
        for (auto it = internal_graph.outcoming_begin(start); it != internal_graph.outcoming_end(start); ++it) {
            auto end = *it;
            if (internal_graph.GetIndegree(end) == 1 and internal_graph.GetOutdegree(start) == 1) {
                result.emplace_back(start, end);
            }
        }
    }
    return result;
}
FilteredReferencePathHelper::FilteredReferencePathHelper(const conj_graph_pack &gp_) : gp_(gp_) {}
FilteredReferencePathHelper::ReferencePaths FilteredReferencePathHelper::GetFilteredReferencePathsFromLength(
        const string &path_to_reference,
        size_t length_threshold) const {
    std::unordered_set<EdgeId> long_edges;
    for (const auto &edge: gp_.g.edges()) {
        if (gp_.g.length(edge) >= length_threshold) {
            long_edges.insert(edge);
        }
    }
    return GetFilteredReferencePathsFromEdges(path_to_reference, long_edges);
}
FilteredReferencePathHelper::ReferencePaths FilteredReferencePathHelper::GetFilteredReferencePathsFromGraph(
        const string &path_to_reference,
        const scaffold_graph::ScaffoldGraph &graph) const {
    //todo implement in the graph
    std::unordered_set<EdgeId> scaffold_edges;
    for (const scaffold_graph::ScaffoldVertex &vertex: graph.vertices()) {
        auto vertex_edges = vertex.GetAllEdges();
        for (const auto &edge: vertex_edges) {
            scaffold_edges.insert(edge);
        }
    }
    return GetFilteredReferencePathsFromEdges(path_to_reference, scaffold_edges);
}
FilteredReferencePathHelper::ReferencePaths FilteredReferencePathHelper::GetFilteredReferencePathsFromEdges(
        const string &path_to_reference,
        const std::unordered_set<EdgeId> &target_edges) const {
    validation::ContigPathBuilder contig_path_builder(gp_);
    auto named_reference_paths = contig_path_builder.GetContigPaths(path_to_reference);
    auto reference_paths = contig_path_builder.StripNames(named_reference_paths);
    DEBUG(reference_paths.size() << " reference paths");
    validation::ContigPathFilter contig_path_filter(target_edges);
    return contig_path_filter.FilterPaths(reference_paths);
}
}
}
}