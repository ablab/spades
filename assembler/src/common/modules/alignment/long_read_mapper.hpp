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
    MappingRange MappingRangeOntoRead_;

    PathWithMappingInfo() = default;
    PathWithMappingInfo(std::vector<EdgeId> && path, MappingRange && range = MappingRange());
};

extern std::mutex PathsWithMappingInfoStorageStorageLock;
extern std::vector<std::unique_ptr<path_extend::BidirectionalPath>> PathsWithMappingInfoStorageStorage;

using PathExtractionF = std::function<std::vector<PathWithMappingInfo> (const MappingPath<EdgeId>&)>;

class LongReadMapper: public SequenceMapperListener {
public:
    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   io::LibraryType lib_type,
                   PathExtractionF path_extractor);

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

        if (lib_type_ == io::LibraryType::TrustedContigs && !paths.empty()) {
            auto path = std::make_unique<path_extend::BidirectionalPath>(g_, paths[0].Path_);
            for (size_t i = 1; i < paths.size(); ++i) {
                auto start_pos = paths[i-1].MappingRangeOntoRead_.initial_range.end_pos;
                start_pos += g_.length(paths[i-1].Path_.back()) - paths[i-1].MappingRangeOntoRead_.mapped_range.end_pos;
                auto end_pos = paths[i].MappingRangeOntoRead_.initial_range.start_pos;
                VERIFY(end_pos >= paths[i].MappingRangeOntoRead_.mapped_range.start_pos);
                end_pos -= paths[i].MappingRangeOntoRead_.mapped_range.start_pos;
                if (end_pos <= start_pos) {
                    std::cout << "start: " << start_pos << " end: " << end_pos << " edgeId: " << paths[i].Path_.front() << '\n';
                }
                std::string gap_seq = (end_pos > start_pos ? r.GetSequenceString().substr(start_pos, end_pos-start_pos) : "");
                VERIFY(g_.k() + end_pos >= start_pos);
                path->PushBack(paths[i].Path_, path_extend::Gap(g_.k() + end_pos - start_pos, std::move(gap_seq)));
            }
            std::lock_guard<std::mutex> guard(PathsWithMappingInfoStorageStorageLock);
            PathsWithMappingInfoStorageStorage.push_back(std::move(path));
        }
        DEBUG("Single read processed");
    }

    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping);

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
    constexpr static double MIN_MAPPED_RATIO = 0.3;
    constexpr static size_t MIN_MAPPED_LENGTH = 100;
public:
    GappedPathExtractor(const Graph& g);

    std::vector<PathWithMappingInfo> operator() (const MappingPath<EdgeId>& mapping) const;

private:

    MappingPath<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const;

    std::vector<PathWithMappingInfo> FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                          MappingPath<EdgeId> &path) const;
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
