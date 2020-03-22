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

struct PathWithMappingInfo { 
    using PathType = std::vector<EdgeId>;
    PathType Path_;
    MappingRange MappingRangeOntoRead_;

    PathWithMappingInfo() = default;
    PathWithMappingInfo(std::vector<EdgeId> && path, MappingRange && range = MappingRange());
};

using PathExtractionF = std::function<std::vector<PathWithMappingInfo> (const MappingPath<EdgeId>&)>;

class LongReadMapper: public SequenceMapperListener {
public:
    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   path_extend::BidirectionalPathStorage& bidirectional_path_storage,
                   io::LibraryType lib_type);

    void StartProcessLibrary(size_t threads_count) override {
        buffer_storages_.reserve(threads_count);
        trusted_path_buffer_storages_.reserve(threads_count);
        for (size_t i = 0; i < threads_count; ++i) {
            buffer_storages_.emplace_back(g_);
            trusted_path_buffer_storages_.emplace_back(g_);
        }
    }

    void StopProcessLibrary() override {
        buffer_storages_.clear();
        trusted_path_buffer_storages_.clear();
    }

    void MergeBuffer(size_t thread_index) override {
        DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index].size());
        storage_.AddStorage(buffer_storages_[thread_index]);
        buffer_storages_[thread_index].Clear();
        std::move(trusted_path_buffer_storages_[thread_index].begin(), trusted_path_buffer_storages_[thread_index].end(), std::back_inserter(trusted_paths_storage_));
        trusted_path_buffer_storages_[thread_index].empty();
        DEBUG("Now size " << storage_.size());
    }

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleRead& r,
                           const MappingPath<EdgeId>& read) override
    {
        ProcessSingleRead(thread_index, read, r);
    }

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleReadSeq&,
                           const MappingPath<EdgeId>& read) override
    {
        ProcessSingleRead(thread_index, read);
    }

    const Graph& g() const {
        return g_;
    }

private:
    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping, const io::SingleRead& r) {
        DEBUG("Processing single read");
        auto paths = path_extractor_(mapping);
        for (const auto& path : paths)
            buffer_storages_[thread_index].AddPath(path.Path_, 1, false);

        auto makeSinglePath = [&]() {
            auto path = std::make_unique<path_extend::BidirectionalPath>(g_, paths[0].Path_);
            for (size_t i = 1; i < paths.size(); ++i) {
                auto start_pos = paths[i-1].MappingRangeOntoRead_.initial_range.end_pos;
                start_pos += g_.length(paths[i-1].Path_.back()) - paths[i-1].MappingRangeOntoRead_.mapped_range.end_pos;

                auto end_pos = paths[i].MappingRangeOntoRead_.initial_range.start_pos;
                VERIFY(end_pos >= paths[i].MappingRangeOntoRead_.mapped_range.start_pos);
                end_pos -= paths[i].MappingRangeOntoRead_.mapped_range.start_pos;

                std::string gap_seq = (end_pos > start_pos ? r.GetSequenceString().substr(start_pos, end_pos-start_pos) : "");

                if ((g_.k() + end_pos < start_pos)) {
                    path.reset();
                    return path;
                }
                path->PushBack(paths[i].Path_, path_extend::Gap(g_.k() + end_pos - start_pos, std::move(gap_seq)));
            }
            return path;
        };

        if (lib_type_ == io::LibraryType::TrustedContigs && !paths.empty()) {
            if (auto path = makeSinglePath()) {
                trusted_path_buffer_storages_[thread_index].push_back(std::move(path));
            } else {
                ERROR("BAD MAPPED PATH WAS DETECTED, DROPPED!");
            }
        }
        DEBUG("Single read processed");
    }

    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping);

    const Graph& g_;
    PathStorage<Graph>& storage_;
    path_extend::BidirectionalPathStorage& trusted_paths_storage_;
    std::vector<PathStorage<Graph>> buffer_storages_;
    std::vector<path_extend::BidirectionalPathStorage> trusted_path_buffer_storages_;
    io::LibraryType lib_type_;
    PathExtractionF path_extractor_;
    DECL_LOGGER("LongReadMapper");
};

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
