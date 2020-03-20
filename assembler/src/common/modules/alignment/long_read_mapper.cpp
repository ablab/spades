//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "long_read_mapper.hpp"

namespace debruijn_graph {

namespace {

class GappedPathExtractor {
public:
    GappedPathExtractor(const Graph& g) noexcept;
    virtual ~GappedPathExtractor() = default;

    std::vector<PathWithMappingInfo> operator() (const MappingPath<EdgeId>& mapping) const;

protected:

    virtual MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const;

    std::vector<PathWithMappingInfo> FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                          MappingPath<EdgeId> &path) const;

    const Graph& g_;
    const MappingPathFixer<Graph> path_fixer_;
    constexpr static double MIN_MAPPED_RATIO = 0.3;
    constexpr static size_t MIN_MAPPED_LENGTH = 100;
};

class GappedPathExtractorForTrustedContigs : public GappedPathExtractor {
public:
    GappedPathExtractorForTrustedContigs(const Graph& g) noexcept;

protected:
    MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const override;
};


std::string GetNum() {
    static std::atomic<unsigned long long> cnt(0);
    return std::to_string(cnt.fetch_add(1));
}

size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path, size_t& mapping_index, MappingRange& range_of_mapped_edge) {
    while(mapping_path[mapping_index].first != edge)
        ++mapping_index;
    size_t start_idx = mapping_index;

    while(mapping_path[mapping_index].first == edge) {
        ++mapping_index;
        if(mapping_index >= mapping_path.size())
            break;
    }
    size_t end_idx = mapping_index;
    size_t total_len = 0;
    for(size_t i = start_idx; i < end_idx; ++i)
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

    if (lib_type == io::LibraryType::PathExtendContigs || lib_type == io::LibraryType::TSLReads 
        || lib_type == io::LibraryType::UntrustedContigs)
    {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractor(g)(mapping);
        };
    }

    return [&] (const MappingPath<EdgeId>& mapping) -> std::vector<PathWithMappingInfo> {
        return {{ReadPathFinder<Graph>(g).FindReadPath(mapping)}};
    };
}

} // namespace

PathWithMappingInfo::PathWithMappingInfo(std::vector<EdgeId> && path, MappingRange && range) 
    : Path_(std::move(path))
    , MappingRangeOntoRead_(range)
{}

LongReadMapper::LongReadMapper(const Graph& g,
                               PathStorage<Graph>& storage,
                               BidirectionalPathStorage& trusted_paths_storage,
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

GappedPathExtractor::GappedPathExtractor(const Graph& g) noexcept
    : g_(g)
    , path_fixer_(g)
{}

GappedPathExtractorForTrustedContigs::GappedPathExtractorForTrustedContigs(const Graph& g) noexcept
    : GappedPathExtractor(g)
{}

std::vector<PathWithMappingInfo> GappedPathExtractor::operator() (const MappingPath<EdgeId>& mapping) const {
    {
        std::ofstream out("mapping_dump" + GetNum()); 
        for (size_t i = 0; i < mapping.size(); ++i)
            out << mapping[i] << '\n';
    }
    auto corrected_path = path_fixer_.DeleteSameEdges(mapping.simple_path());
    auto filtered_path = FilterBadMappings(corrected_path, mapping);
    return FindReadPathWithGaps(mapping, filtered_path);
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
        if (mapping_size > MIN_MAPPED_LENGTH || 
            math::gr((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO))
        {
            new_corrected_path.push_back(edge, range_of_mapped_edge);
        }
    }
    return new_corrected_path;
}

MappingPath<EdgeId> GappedPathExtractorForTrustedContigs::FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                                            const MappingPath<EdgeId> &mapping_path) const
{
    MappingPath<EdgeId> new_corrected_path;
    size_t mapping_index = 0;
    for (auto edge : corrected_path) {
        MappingRange range_of_mapped_edge;
        size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index, range_of_mapped_edge);
        size_t edge_len =  g_.length(edge);
        if (math::ls((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO + 1) &&
            math::gr((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO))
        {
            if (!new_corrected_path.empty() && new_corrected_path.back().first == edge && new_corrected_path.back().second.mapped_range.end_pos <= range_of_mapped_edge.mapped_range.start_pos) {
                range_of_mapped_edge.mapped_range.start_pos = new_corrected_path.back().second.mapped_range.start_pos;
                new_corrected_path.pop_back();
            }
            new_corrected_path.push_back(edge, range_of_mapped_edge);
        }
    }
    return new_corrected_path;
}

std::vector<PathWithMappingInfo> GappedPathExtractor::FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                                           MappingPath<EdgeId> &path) const
{
    std::vector<PathWithMappingInfo> result;
    if (mapping_path.empty()) {
        TRACE("read unmapped");
        return result;
    }
    std::ofstream out("path_dump" + GetNum()); 
    out << std::setw(8) << path[0].first << ' ' << path[0].second << '\n';
    PathWithMappingInfo tmp_path;
    auto setStartPos = [&tmp_path] (auto const & path) {
        tmp_path.MappingRangeOntoRead_.initial_range.start_pos = path.second.initial_range.start_pos;
        tmp_path.MappingRangeOntoRead_.mapped_range.start_pos = path.second.mapped_range.start_pos;
    };
    auto setEndPos = [&tmp_path] (auto const & path) {
        tmp_path.MappingRangeOntoRead_.initial_range.end_pos = path.second.initial_range.end_pos;
        tmp_path.MappingRangeOntoRead_.mapped_range.end_pos = path.second.mapped_range.end_pos;
    };

    setStartPos(path[0]);
    tmp_path.Path_.push_back(path[0].first);
    for (size_t i = 1; i < path.size(); ++i) {
        auto left_vertex = g_.EdgeEnd(path[i - 1].first);
        auto right_vertex = g_.EdgeStart(path[i].first);
        if (left_vertex != right_vertex) {
            auto closure = path_fixer_.TryCloseGap(left_vertex, right_vertex);
            if (!closure.empty()) {
                tmp_path.Path_.insert(tmp_path.Path_.end(), closure.begin(), closure.end());
                out << "gap is closed!\n";
            } else {
                setEndPos(path[i-1]);
                out << "gap is not closed, length of path: " << tmp_path.MappingRangeOntoRead_.initial_range.end_pos - tmp_path.MappingRangeOntoRead_.initial_range.start_pos << '\n';
                result.push_back(std::move(tmp_path));
                setStartPos(path[i]);
            }
        }
        // VERIFY(tmp_path.Path_.empty() || bool(path[i-1].second.end_pos == path[i].second.start_pos));
        if (!(tmp_path.Path_.empty() || bool(path[i-1].second.initial_range.end_pos == path[i].second.initial_range.start_pos)))
            out << "OH NO!\n";
        tmp_path.Path_.push_back(path[i].first);
        out << std::setw(8) << path[i].first << ' ' << path[i].second << '\n';
    }
    setEndPos(path.back());
    out << "length of path: " << tmp_path.MappingRangeOntoRead_.initial_range.end_pos - tmp_path.MappingRangeOntoRead_.initial_range.start_pos << '\n';
    result.push_back(std::move(tmp_path));

    for (int i = 1; i < result.size(); ++i)
        VERIFY(result[i-1].MappingRangeOntoRead_.initial_range.end_pos < result[i].MappingRangeOntoRead_.initial_range.start_pos);

    std::cout << "================DONE====================" << std::endl; 
    return result;
}

} // namespace debruijn_graph

