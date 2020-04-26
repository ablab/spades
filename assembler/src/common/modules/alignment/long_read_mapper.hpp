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

#include "pipeline/library_fwd.hpp"

namespace debruijn_graph {

typedef std::function<std::vector<std::vector<EdgeId>> (const MappingPath<EdgeId>&)> PathExtractionF;

class LongReadMapper: public SequenceMapperListener {
public:
    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   PathExtractionF path_extractor);

    void StartProcessLibrary(size_t threads_count) override;

    void StopProcessLibrary() override;

    void MergeBuffer(size_t thread_index) override;

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleRead&,
                           const MappingPath<EdgeId>& read) override;

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleReadSeq&,
                           const MappingPath<EdgeId>& read) override;

    const Graph& g() const {
        return g_;
    }

private:
    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping);

    const Graph& g_;
    PathStorage<Graph>& storage_;
    std::vector<PathStorage<Graph>> buffer_storages_;
    PathExtractionF path_extractor_;
    DECL_LOGGER("LongReadMapper");
};

class GappedPathExtractor {
    const Graph& g_;
    const MappingPathFixer<Graph> path_fixer_;
    const double MIN_MAPPED_RATIO = 0.3;
    const size_t MIN_MAPPED_LENGTH = 100;
public:
    GappedPathExtractor(const Graph& g);

    std::vector<std::vector<EdgeId>> operator() (const MappingPath<EdgeId>& mapping) const;

private:

    size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path, size_t& mapping_index) const;

    std::vector<EdgeId> FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                          const MappingPath<EdgeId> &mapping_path) const;

    std::vector<std::vector<EdgeId>> FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                          std::vector<EdgeId> &corrected_path) const;

    std::vector<std::vector<EdgeId>> SplitUnfixedPoints(std::vector<EdgeId>& path) const;
};

PathExtractionF ChooseProperReadPathExtractor(const Graph& g, io::LibraryType lib_type);

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
