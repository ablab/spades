//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "transition_extractor.hpp"

#include "io/reads/file_reader.hpp"
#include "io/reads/rc_reader_wrapper.hpp"
#include "alignment/sequence_mapper.hpp"
#include "alignment/long_read_mapper.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

//fixme duplication with long_read_mapper

size_t GappedPathExtractor::CountMappedEdgeSize(EdgeId edge,
                                                const MappingPath<EdgeId> &mapping_path,
                                                size_t &mapping_index,
                                                MappingRange &range_of_mapped_edge) const
{
    while (mapping_path[mapping_index].first != edge)
        ++mapping_index;
    size_t start_idx = mapping_index;

    while (mapping_path[mapping_index].first == edge) {
        ++mapping_index;
        if (mapping_index >= mapping_path.size())
            break;
    }
    size_t end_idx = mapping_index;
    size_t total_len = 0;
    for (size_t i = start_idx; i < end_idx; ++i)
        total_len += mapping_path[i].second.initial_range.size();

    range_of_mapped_edge.initial_range.start_pos = mapping_path[start_idx].second.initial_range.start_pos;
    range_of_mapped_edge.mapped_range.start_pos = mapping_path[start_idx].second.mapped_range.start_pos;
    range_of_mapped_edge.initial_range.end_pos = mapping_path[end_idx-1].second.initial_range.end_pos;
    range_of_mapped_edge.mapped_range.end_pos = mapping_path[end_idx-1].second.mapped_range.end_pos;

    return total_len;
}

GappedPathExtractor::GappedPathExtractor(const Graph& g) noexcept
    : g_(g)
    , path_fixer_(g)
{}

std::vector<PathWithMappingInfo> GappedPathExtractor::operator() (const MappingPath<EdgeId>& mapping) const {
    auto corrected_path = path_fixer_.DeleteSameEdges(mapping.simple_path());
    auto filtered_path = FilterBadMappings(corrected_path, mapping);
    return FindReadPathWithGaps(filtered_path);
}

MappingPath<EdgeId> GappedPathExtractor::FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                           const MappingPath<EdgeId> &mapping_path) const
{
    MappingPath<EdgeId> new_corrected_path;
    size_t mapping_index = 0;
    for (auto edge : corrected_path) {
        MappingRange range_of_mapped_edge;
        size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index, range_of_mapped_edge);
        size_t edge_len =  g_.length(edge);
        if (mapping_size > MIN_MAPPED_LENGTH || math::gr((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO))
            new_corrected_path.push_back(edge, range_of_mapped_edge);
    }
    return new_corrected_path;
}

std::vector<PathWithMappingInfo> GappedPathExtractor::FindReadPathWithGaps(MappingPath<EdgeId> &path) const {
    std::vector<PathWithMappingInfo> result;
    if (path.empty()) {
        TRACE("read unmapped");
        return result;
    }
    PathWithMappingInfo tmp_path;
    auto SetStartPos = [&tmp_path] (auto const & path) {
      tmp_path.MappingRangeOntoRead_.initial_range.start_pos = path.second.initial_range.start_pos;
      tmp_path.MappingRangeOntoRead_.mapped_range.start_pos = path.second.mapped_range.start_pos;
    };
    auto SetEndPos = [&tmp_path] (auto const & path) {
      tmp_path.MappingRangeOntoRead_.initial_range.end_pos = path.second.initial_range.end_pos;
      tmp_path.MappingRangeOntoRead_.mapped_range.end_pos = path.second.mapped_range.end_pos;
    };

    SetStartPos(path[0]);
    tmp_path.Path_.push_back(path[0].first);
    for (size_t i = 1; i < path.size(); ++i) {
        auto left_vertex = g_.EdgeEnd(path[i - 1].first);
        auto right_vertex = g_.EdgeStart(path[i].first);
        if (left_vertex != right_vertex) {
            auto closure = path_fixer_.TryCloseGap(left_vertex, right_vertex);
            if (!closure.empty()) {
                tmp_path.Path_.insert(tmp_path.Path_.end(), closure.begin(), closure.end());
            } else {
                SetEndPos(path[i-1]);
                result.push_back(std::move(tmp_path));
                SetStartPos(path[i]);
            }
        }
        tmp_path.Path_.push_back(path[i].first);
    }
    SetEndPos(path.back());
    result.push_back(std::move(tmp_path));
    return result;
}

//duplication end

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
    io::FileReadStream contig_stream(contig_path);
    auto rc_contig_stream_ptr = std::make_shared<io::RCWrapper<io::SingleRead>>(std::move(contig_stream));

    std::vector<io::SingleRead> contigs;
    const size_t contig_length_threshold = path_length_threshold_;
    while (!rc_contig_stream_ptr->eof()) {
        io::SingleRead contig;
        (*rc_contig_stream_ptr) >> contig;
        auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, debruijn_graph::EdgeIndex<Graph>>>
            (g_, index_, kmer_mapper_);
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
    GappedPathExtractor gapped_path_extractor(g_);
    for (const auto &named_path: contig_paths) {
        auto path = named_path.mapping_path;
        auto fixed_paths = gapped_path_extractor(path);
        std::vector<EdgeWithMapping> long_path;
        size_t prefix_len = 0;
        for (const auto &fixed_path: fixed_paths) {
            for (const auto &edge: fixed_path.Path_) {
                Range mapping(prefix_len, prefix_len + g_.length(edge));
                EdgeWithMapping ewm(edge, mapping);
                long_path.push_back(ewm);
                prefix_len += g_.length(edge);
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
            if (pred_(ewm.edge_)) {
                filtered_path.push_back(ewm);
            }
        }
        if (not filtered_path.empty()) {
            filtered_paths.push_back(filtered_path);
        }
    }
    auto merged_paths = MergeSameEdges(filtered_paths);

//    auto result = RemoveRepeats(merged_paths);

    return merged_paths;
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

ContigTransitionStorage StrictTransitionStorageBuilder::BuildStorage(const UniqueReferencePaths &contig_paths) const {
    ContigTransitionStorage storage(contig_paths.unique_storage_);
    for (const auto &path: contig_paths.paths_) {
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
ContigTransitionStorage ReverseTransitionStorageBuilder::BuildStorage(const UniqueReferencePaths &contig_paths) const {
    ContigTransitionStorage storage(contig_paths.unique_storage_);
    for (const auto &path: contig_paths.paths_) {
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
ContigTransitionStorage ConjugateTransitionStorageBuilder::BuildStorage(const UniqueReferencePaths &contig_paths) const {
    ContigTransitionStorage storage(contig_paths.unique_storage_);
    for (const auto &path: contig_paths.paths_) {
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
ContigTransitionStorage GeneralTransitionStorageBuilder::BuildStorage(const UniqueReferencePaths &contig_paths) const {
    ContigTransitionStorage storage(contig_paths.unique_storage_);
    for (const auto &path: contig_paths.paths_) {
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
std::vector<ClusterTransitionExtractor::Transition> ClusterTransitionExtractor::ExtractAllTransitionsFromNonPathCluster(
    const cluster_storage::Cluster &cluster) {
    const auto &internal_graph = cluster.GetInternalGraph();
    std::vector<Transition> result;
    for (const auto &vertex: internal_graph) {
        for (const auto &other: internal_graph.OutNeighbours(vertex)) {
            result.emplace_back(vertex, other);
        }
    }
    return result;
}
std::vector<ClusterTransitionExtractor::Transition> ClusterTransitionExtractor::ExtractGoodTransitionsFromNonPathCluster(
    const cluster_storage::Cluster &cluster) {
    const auto &internal_graph = cluster.GetInternalGraph();
    std::vector<Transition> result;
    for (const auto &start: internal_graph) {
        for (const auto &end: internal_graph.OutNeighbours(start)) {
            if (internal_graph.GetIndegree(end) == 1 and internal_graph.GetOutdegree(start) == 1) {
                result.emplace_back(start, end);
            }
        }
    }
    return result;
}
UniqueReferencePaths FilteredReferencePathHelper::GetFilteredReferencePathsFromUnique(
    const string &path_to_reference,
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage) const {
    validation::ContigPathBuilder contig_path_builder(g_, index_, kmer_mapper_, unique_storage.min_length());
    auto named_reference_paths = contig_path_builder.GetContigPaths(path_to_reference);
    auto reference_paths = contig_path_builder.StripNames(named_reference_paths);
    DEBUG(reference_paths.size() << " reference paths");
    auto is_unique = [&unique_storage](const EdgeId &edge) {
      return unique_storage.IsUnique(edge);
    };
    validation::ContigPathFilter contig_path_filter(is_unique);
    auto paths = contig_path_filter.FilterPaths(reference_paths);
    UniqueReferencePaths result(std::move(paths), unique_storage);
    return result;
}
}
}
}