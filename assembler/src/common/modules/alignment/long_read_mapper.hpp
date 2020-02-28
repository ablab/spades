//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "long_read_storage.hpp"
#include "sequence_mapper_notifier.hpp"

#include <cassert>
#include <mutex>
#include <iomanip>

namespace debruijn_graph {

//todo: extend interface of Path<EdgeId> and replace to its
struct PathWithMappingInfo { 
    using PathType = std::vector<EdgeId>;
    PathType Path_;
    Range MappingRangeOntoRead_;
    uint64_t Id_ = -1;

    PathWithMappingInfo() = default;
    PathWithMappingInfo(std::vector<EdgeId> && path, Range && range = Range());
};

template<class ReadType>
struct PathsWithMappingInfoStorage {
    ReadType Read_;
    std::vector<PathWithMappingInfo> Paths_;
};

extern std::mutex PathsWithMappingInfoStorageStorageLock;
extern std::vector<PathsWithMappingInfoStorage<io::SingleRead>> PathsWithMappingInfoStorageStorage;

using PathExtractionF = std::function<std::vector<PathWithMappingInfo> (const MappingPath<EdgeId>&)>;

class LongReadMapper: public SequenceMapperListener {
public:
    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   io::LibraryType lib_type,
                   PathExtractionF path_extractor)
            : g_(g),
              storage_(storage),
              lib_type_(lib_type),
              path_extractor_(path_extractor) {
    }

    void StartProcessLibrary(size_t threads_count) override;

    void StopProcessLibrary() override;

    void MergeBuffer(size_t thread_index) override;

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleRead& r,
                           const MappingPath<EdgeId>& read) override
    {
        ProcessSingleRead(thread_index, read, r);
    }

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleReadSeq&,
                           const MappingPath<EdgeId>& read) override;

    const Graph& g() const {
        return g_;
    }

private:
    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping, const io::SingleRead& r) {
        DEBUG("Processing single read");
        auto paths = path_extractor_(mapping);
        for (const auto& path : paths)
            buffer_storages_[thread_index].AddPath(path.Path_, 1, false);

        if (lib_type_ == io::LibraryType::TrustedContigs) {
            std::lock_guard<std::mutex> guard(PathsWithMappingInfoStorageStorageLock);
            PathsWithMappingInfoStorageStorage.push_back({r, std::move(paths)});
        }
        DEBUG("Single read processed");
    }

    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping) {
        DEBUG("Processing read");
        for (const auto& path : path_extractor_(mapping)) {
            buffer_storages_[thread_index].AddPath(path.Path_, 1, false);
        }
        DEBUG("Read processed");
    }

    const Graph& g_;
    PathStorage<Graph>& storage_;
    std::vector<PathStorage<Graph>> buffer_storages_;
    io::LibraryType lib_type_;
    PathExtractionF path_extractor_;
    DECL_LOGGER("LongReadMapper");
};

class GappedPathExtractor {
    const Graph& g_;
    const MappingPathFixer<Graph> path_fixer_;
    const double MIN_MAPPED_RATIO = 0.3;
    const size_t MIN_MAPPED_LENGTH = 100;
    using MappedEdge = std::pair<EdgeId, Range>;
public:
    GappedPathExtractor(const Graph& g);

    std::vector<PathWithMappingInfo> operator() (const MappingPath<EdgeId>& mapping) const {
        auto corrected_path = path_fixer_.DeleteSameEdges(mapping.simple_path());
        auto filtered_path = FilterBadMappings(corrected_path, mapping);
        return FindReadPathWithGaps(mapping, filtered_path);
    }

private:

    size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path, size_t& mapping_index, Range& range_of_mapped_edge) const {
        while(mapping_path[mapping_index].first != edge) {
            mapping_index++;
        }
        size_t start_idx = mapping_index;

        while(mapping_path[mapping_index].first == edge) {
            mapping_index++;
            if(mapping_index >= mapping_path.size()) {
                break;
            }
        }
        size_t end_idx = mapping_index;
        size_t total_len = 0;
        for(size_t i = start_idx; i < end_idx; ++i) {
            total_len += mapping_path[i].second.initial_range.size();
        }

        range_of_mapped_edge.start_pos = mapping_path[start_idx].second.initial_range.start_pos;
        range_of_mapped_edge.end_pos = mapping_path[end_idx-1].second.initial_range.end_pos;

        return total_len;
    }

    std::vector<MappedEdge> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const {
        std::vector<MappedEdge> new_corrected_path;
        size_t mapping_index = 0;
        for (auto edge : corrected_path) {
            Range range_of_mapped_edge;
            size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index, range_of_mapped_edge);
            size_t edge_len =  g_.length(edge);
            if (mapping_size > MIN_MAPPED_LENGTH || 
                    math::gr((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO)) {
                new_corrected_path.emplace_back(edge, range_of_mapped_edge);
            }
        }
        return new_corrected_path;
    }

    std::vector<PathWithMappingInfo> FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                          std::vector<MappedEdge> &path) const
    {
        std::vector<PathWithMappingInfo> result;
        if (mapping_path.empty()) {
            TRACE("read unmapped");
            return result;
        }
        std::ofstream out("path_dump"); 
        out << std::setw(8) << path[0].first << ' ' << path[0].second << '\n';
        PathWithMappingInfo tmp_path;
        tmp_path.Path_.push_back(path[0].first);
        tmp_path.MappingRangeOntoRead_.start_pos = path[0].second.start_pos;
        for (size_t i = 1; i < path.size(); ++i) {
            auto left_vertex = g_.EdgeEnd(path[i - 1].first);
            auto right_vertex = g_.EdgeStart(path[i].first);
            if (left_vertex != right_vertex) {
                auto closure = path_fixer_.TryCloseGap(left_vertex, right_vertex);
                if (!closure.empty()) {
                    tmp_path.Path_.insert(tmp_path.Path_.end(), closure.begin(), closure.end());
                    out << "gap is closed!\n";
                } else {
                    tmp_path.MappingRangeOntoRead_.end_pos = path[i-1].second.end_pos;
                    result.push_back(std::move(tmp_path));
                    tmp_path.MappingRangeOntoRead_.start_pos = path[i].second.start_pos;
                    out << "gap is not closed\n";
                }
            }
            VERIFY(tmp_path.Path_.empty() || bool(path[i-1].second.end_pos == path[i].second.start_pos));
            tmp_path.Path_.push_back(path[i].first);
            out << std::setw(8) << path[i].first << ' ' << path[i].second << '\n';
        }
        tmp_path.MappingRangeOntoRead_.end_pos = path.back().second.end_pos;
        result.push_back(std::move(tmp_path));

        for (int i = 1; i < result.size(); ++i)
            VERIFY(result[i-1].MappingRangeOntoRead_.end_pos < result[i].MappingRangeOntoRead_.start_pos);

        std::cout << "================DONE====================" << std::endl; 
        return result;
    }

};

inline PathExtractionF ChooseProperReadPathExtractor(const Graph& g, io::LibraryType lib_type) {
    if (lib_type == io::LibraryType::PathExtendContigs || lib_type == io::LibraryType::TSLReads
        || lib_type == io::LibraryType::TrustedContigs || lib_type == io::LibraryType::UntrustedContigs) {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractor(g)(mapping);
        };
    } else {
        return [&] (const MappingPath<EdgeId>& mapping) -> std::vector<PathWithMappingInfo> {
            return {{ReadPathFinder<Graph>(g).FindReadPath(mapping)}};
        };
    }
}

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
