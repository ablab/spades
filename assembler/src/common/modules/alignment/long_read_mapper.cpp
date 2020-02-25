//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "long_read_mapper.hpp"

namespace debruijn_graph {

LongReadMapper::LongReadMapper(const Graph& g,
                               PathStorage<Graph>& storage,
                               PathExtractionF path_extractor)
    : g_(g)
    , storage_(storage)
    , path_extractor_(path_extractor) 
{}

void LongReadMapper::StartProcessLibrary(size_t threads_count) {
    for (size_t i = 0; i < threads_count; ++i)
        buffer_storages_.emplace_back(g_);
}

void LongReadMapper::StopProcessLibrary() {
    buffer_storages_.clear();
}

void LongReadMapper::MergeBuffer(size_t thread_index) {
    DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index].size());
    storage_.AddStorage(buffer_storages_[thread_index]);
    buffer_storages_[thread_index].Clear();
    DEBUG("Now size " << storage_.size());
}

void LongReadMapper::ProcessSingleRead(size_t thread_index,
                                       const io::SingleRead&,
                                       const MappingPath<EdgeId>& read) 
{
    ProcessSingleRead(thread_index, read);
}

void LongReadMapper::ProcessSingleRead(size_t thread_index,
                                       const io::SingleReadSeq&,
                                       const MappingPath<EdgeId>& read)
{
    ProcessSingleRead(thread_index, read);
}

void LongReadMapper::ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping) {
    DEBUG("Processing read");
    for (const auto& path : path_extractor_(mapping)) {
        buffer_storages_[thread_index].AddPath(path, 1, false);
    }
    DEBUG("Read processed");
}

GappedPathExtractor::GappedPathExtractor(const Graph& g)
    : g_(g)
    , path_fixer_(g) 
{}

std::vector<std::vector<EdgeId>> GappedPathExtractor::operator() (const MappingPath<EdgeId>& mapping) const {
    auto corrected_path = path_fixer_.DeleteSameEdges(mapping.simple_path());
    corrected_path = FilterBadMappings(corrected_path, mapping);
    return FindReadPathWithGaps(mapping, corrected_path);
}

size_t GappedPathExtractor::CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path, size_t& mapping_index) const {
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

    return total_len;
}

std::vector<EdgeId> GappedPathExtractor::FilterBadMappings(const std::vector<EdgeId> &corrected_path,
                                                           const MappingPath<EdgeId> &mapping_path) const
{
    std::vector<EdgeId> new_corrected_path;
    size_t mapping_index = 0;
    for (auto edge : corrected_path) {
        size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index);
        size_t edge_len =  g_.length(edge);
        //VERIFY(edge_len >= mapping_size);
        if (mapping_size > MIN_MAPPED_LENGTH || 
                math::gr((double) mapping_size / (double) edge_len, MIN_MAPPED_RATIO)) {
            new_corrected_path.push_back(edge);
        }
    }
    return new_corrected_path;
}

std::vector<std::vector<EdgeId>> GappedPathExtractor::FindReadPathWithGaps(const MappingPath<EdgeId> &mapping_path,
                                                                           std::vector<EdgeId> &corrected_path) const 
{
    if (mapping_path.size() == 0) {
        TRACE("read unmapped");
        return {};
    }
    auto fixed_path = path_fixer_.TryFixPath(corrected_path);
    return SplitUnfixedPoints(fixed_path);
}

std::vector<std::vector<EdgeId>> GappedPathExtractor::SplitUnfixedPoints(std::vector<EdgeId>& path) const {
    std::vector<std::vector<EdgeId>> result;
    size_t prev_start = 0;
    for (size_t i = 1; i < path.size(); ++i) {
        if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i])) {
                result.emplace_back(path.begin() + prev_start, path.begin() + i);
                prev_start = i;
        }
    }
    result.emplace_back(path.begin() + prev_start, path.end());
    return result;
}

PathExtractionF ChooseProperReadPathExtractor(const Graph& g, io::LibraryType lib_type) {
    if (lib_type == io::LibraryType::PathExtendContigs || lib_type == io::LibraryType::TSLReads
        || lib_type == io::LibraryType::TrustedContigs || lib_type == io::LibraryType::UntrustedContigs) {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractor(g)(mapping);
        };
    } else {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return std::vector<std::vector<EdgeId>>{ReadPathFinder<Graph>(g).FindReadPath(mapping)};
        };
    }
}

} // namespace debruijn_graph