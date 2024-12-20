//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "long_read_storage.hpp"
#include "sequence_mapper_notifier.hpp"

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "library/library_fwd.hpp"

namespace debruijn_graph {
// FIXME: get rid of this
using omnigraph::MappingPath;
using omnigraph::MappingRange;

struct PathWithMappingInfo {
    using PathType = std::vector<EdgeId>;
    PathType Path_;
    omnigraph::MappingRange MappingRangeOntoRead_;

    PathWithMappingInfo() = default;
    PathWithMappingInfo(std::vector<EdgeId> && path, omnigraph::MappingRange && range = omnigraph::MappingRange());
};

using PathExtractionF = std::function<std::vector<PathWithMappingInfo> (const omnigraph::MappingPath<EdgeId>&)>;

class LongReadMapper: public SequenceMapperListener {
  public:
    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   path_extend::GappedPathStorage& bidirectional_path_storage,
                   io::LibraryType lib_type);

    void StartProcessLibrary(size_t threads_count) override;
    void StopProcessLibrary() override;

    void MergeBuffer(size_t thread_index) override;

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleRead& r,
                           const omnigraph::MappingPath<EdgeId>& read) override;

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleReadSeq&,
                           const omnigraph::MappingPath<EdgeId>& read) override;

    const Graph& g() const noexcept {
        return g_;
    }

    void Serialize(std::ostream &os) const override {
        storage_.BinWrite(os);
    }

    void Deserialize(std::istream &is) override {
        storage_.BinRead(is);
    }

    void MergeFromStream(std::istream &is) override {
        PathStorage<Graph> remote(g_);
        remote.BinRead(is);
        storage_.AddStorage(remote);
    }

private:

    void ProcessSingleRead(size_t thread_index, const omnigraph::MappingPath<EdgeId>& mapping, const io::SingleRead& r);

    void ProcessSingleRead(size_t thread_index, const omnigraph::MappingPath<EdgeId>& mapping);

    const Graph& g_;
    PathStorage<Graph>& storage_;
    path_extend::GappedPathStorage& trusted_paths_storage_;
    std::vector<PathStorage<Graph>> buffer_storages_;
    std::vector<path_extend::GappedPathStorage> trusted_path_buffer_storages_;
    io::LibraryType lib_type_;
    PathExtractionF path_extractor_;
    DECL_LOGGER("LongReadMapper");
};

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
