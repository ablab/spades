//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "long_read_mapper.hpp"

namespace debruijn_graph {

namespace {

/// @brief Make a set of paths with mapping info
class GappedPathExtractor {
public:
    GappedPathExtractor(const Graph& g) noexcept;
    virtual ~GappedPathExtractor() = default;

    std::vector<PathWithMappingInfo> operator() (const MappingPath<EdgeId>& mapping) const;

protected:

    virtual MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                  const MappingPath<EdgeId> &mapping_path) const;

    std::vector<PathWithMappingInfo> FindReadPathWithGaps(MappingPath<EdgeId> &path) const;

    const Graph& g_;
    const MappingPathFixer<Graph> path_fixer_;
    const static double MIN_MAPPED_RATIO;
    constexpr static size_t MIN_MAPPED_LENGTH = 100;
};

const double GappedPathExtractor::MIN_MAPPED_RATIO = 0.3;

class GappedPathExtractorForTrustedContigs : public GappedPathExtractor {
public:
    GappedPathExtractorForTrustedContigs(const Graph& g) noexcept;

protected:
    MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const override;
};

/// finds the first subsequence in 'mapping_path' starting from 'mapping_index' witch contains only edges equal to 'edge'.
/// @returns the sum of the edge lengths of this subsequence.
/// sets mapping range of this subsequence to the 'range_of_mapped_edge'.
size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path,
                           size_t& mapping_index, MappingRange& range_of_mapped_edge)
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

PathExtractionF ChooseProperReadPathExtractor(const Graph& g, io::LibraryType lib_type) {
    if (lib_type == io::LibraryType::TrustedContigs) {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractorForTrustedContigs(g)(mapping);
        };
    }

    if (lib_type == io::LibraryType::PathExtendContigs || 
        lib_type == io::LibraryType::TSLReads ||
        lib_type == io::LibraryType::UntrustedContigs) {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractor(g)(mapping);
        };
    }

    return [&] (const MappingPath<EdgeId>& mapping) -> std::vector<PathWithMappingInfo> {
        return {{ReadPathFinder<Graph>(g).FindReadPath(mapping)}};
    };
}

template<class T>
inline bool InBounds(T const & min_value, T const & middle_value, T const & max_value) {
    return math::ls(min_value, middle_value) && math::ls(middle_value, max_value);
}

/// merges 'paths' into a single path, filling gaps with contig substrings
std::vector<path_extend::SimpleBidirectionalPath> MergePaths (const std::vector<PathWithMappingInfo> & paths, const Graph & g_, const io::SingleRead& r) {
    std::vector<path_extend::SimpleBidirectionalPath> merged_paths;
    merged_paths.push_back(path_extend::SimpleBidirectionalPath(paths[0].Path_));
    for (size_t i = 1; i < paths.size(); ++i) {
        auto start_pos = paths[i-1].MappingRangeOntoRead_.initial_range.end_pos;
        start_pos += g_.length(paths[i-1].Path_.back()) - paths[i-1].MappingRangeOntoRead_.mapped_range.end_pos;

        size_t end_pos = paths[i].MappingRangeOntoRead_.initial_range.start_pos;
        end_pos -= paths[i].MappingRangeOntoRead_.mapped_range.start_pos;

        std::string gap_seq = (end_pos > start_pos ? r.GetSequenceString().substr(start_pos, end_pos-start_pos) : "");

        if ((g_.k() + end_pos < start_pos)) {
            merged_paths.push_back(path_extend::SimpleBidirectionalPath(paths[i].Path_));
        } else {
            auto k = static_cast<int>(g_.k() + end_pos - start_pos);
            merged_paths.back().PushBack(paths[i].Path_, path_extend::Gap(std::move(gap_seq), k));
        }
    }
    return merged_paths;
};


} // namespace

PathWithMappingInfo::PathWithMappingInfo(std::vector<EdgeId> && path, MappingRange && range) 
    : Path_(std::move(path))
    , MappingRangeOntoRead_(range)
{}

LongReadMapper::LongReadMapper(const Graph& g,
                               PathStorage<Graph>& storage,
                               path_extend::GappedPathStorage& trusted_paths_storage,
                               io::LibraryType lib_type)
    : g_(g)
    , storage_(storage)
    , trusted_paths_storage_(trusted_paths_storage)
    , lib_type_(lib_type)
    , path_extractor_(ChooseProperReadPathExtractor(g, lib_type))
{}

void LongReadMapper::ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping) {
    DEBUG("Processing read");
    for (const auto& path : path_extractor_(mapping))
        buffer_storages_[thread_index].AddPath(path.Path_, 1, false);
    DEBUG("Read processed");
}

void LongReadMapper::StartProcessLibrary(size_t threads_count) {
    buffer_storages_.reserve(threads_count);
    trusted_path_buffer_storages_.reserve(threads_count);
    for (size_t i = 0; i < threads_count; ++i) {
        buffer_storages_.emplace_back(g_);
        trusted_path_buffer_storages_.emplace_back();
    }
}

void LongReadMapper::StopProcessLibrary() {
    buffer_storages_.clear();
    trusted_path_buffer_storages_.clear();
}

void LongReadMapper::MergeBuffer(size_t thread_index) {
    DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index].size());
    storage_.AddStorage(buffer_storages_[thread_index]);
    buffer_storages_[thread_index].Clear();
    std::move(trusted_path_buffer_storages_[thread_index].begin(), trusted_path_buffer_storages_[thread_index].end(), std::back_inserter(trusted_paths_storage_));
    trusted_path_buffer_storages_[thread_index].clear();
    DEBUG("Now size " << storage_.size());
}

void LongReadMapper::ProcessSingleRead(size_t thread_index,
                                       const io::SingleRead& r,
                                       const MappingPath<EdgeId>& read)
{
    ProcessSingleRead(thread_index, read, r);
}

void LongReadMapper::ProcessSingleRead(size_t thread_index,
                                       const io::SingleReadSeq&,
                                       const MappingPath<EdgeId>& read)
{
    ProcessSingleRead(thread_index, read);
}

void LongReadMapper::ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping, const io::SingleRead& r) {
    DEBUG("Processing single read");
    auto paths = path_extractor_(mapping);
    for (const auto& path : paths)
        buffer_storages_[thread_index].AddPath(path.Path_, 1, false);

    if (lib_type_ == io::LibraryType::TrustedContigs && !paths.empty()) {
        auto gapped_paths = MergePaths(paths, g_, r);
        std::move(gapped_paths.begin(), gapped_paths.end(), std::back_inserter(trusted_path_buffer_storages_[thread_index]));
    }
    DEBUG("Single read processed");
}

GappedPathExtractor::GappedPathExtractor(const Graph& g) noexcept
    : g_(g)
    , path_fixer_(g)
{}

GappedPathExtractorForTrustedContigs::GappedPathExtractorForTrustedContigs(const Graph& g) noexcept
    : GappedPathExtractor(g)
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

MappingPath<EdgeId> GappedPathExtractorForTrustedContigs::FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                                            const MappingPath<EdgeId> &mapping_path) const
{
    MappingPath<EdgeId> new_corrected_path;
    MappingRange range_of_mapped_edge;
    auto is_the_previous_edge_extension = [&new_corrected_path, &range_of_mapped_edge](EdgeId edge) {
        return !new_corrected_path.empty() &&
                new_corrected_path.back().first == edge &&
                new_corrected_path.back().second.mapped_range.end_pos <= range_of_mapped_edge.mapped_range.start_pos;
    };
    size_t mapping_index = 0;
    for (auto edge : corrected_path) {
        size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index, range_of_mapped_edge);
        size_t edge_len =  g_.length(edge);
        if (InBounds(MIN_MAPPED_RATIO, (double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO + 1)) {
            if (is_the_previous_edge_extension(edge)) {
                range_of_mapped_edge.mapped_range.start_pos = new_corrected_path.back().second.mapped_range.start_pos;
                new_corrected_path.pop_back();
            }
            new_corrected_path.push_back(edge, range_of_mapped_edge);
        }
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

} // namespace debruijn_graph

