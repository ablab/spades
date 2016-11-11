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

namespace debruijn_graph {

class LongReadMapper: public SequenceMapperListener {
public:
    typedef vector<vector<EdgeId>> PathsT;
    typedef MappingPath<EdgeId> MappingT;
    typedef std::function<PathsT (const MappingT&)> PathExtractionF;

    LongReadMapper(const Graph& g,
                   PathStorage<Graph>& storage,
                   PathExtractionF path_extractor)
            : g_(g),
              storage_(storage),
              path_extractor_(path_extractor) {
    }

    void StartProcessLibrary(size_t threads_count) override {
        for (size_t i = 0; i < threads_count; ++i)
            buffer_storages_.emplace_back(g_);
    }

    void StopProcessLibrary() override {
        buffer_storages_.clear();
    }

    void MergeBuffer(size_t thread_index) override {
        DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index].size());
        storage_.AddStorage(buffer_storages_[thread_index]);
        buffer_storages_[thread_index].Clear();
        DEBUG("Now size " << storage_.size());
    }

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleRead&,
                           const MappingPath<EdgeId>& read) override {
        ProcessSingleRead(thread_index, read);
    }

    void ProcessSingleRead(size_t thread_index,
                           const io::SingleReadSeq&,
                           const MappingPath<EdgeId>& read) override {
        ProcessSingleRead(thread_index, read);
    }

    const Graph& g() const {
        return g_;
    }

private:
    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& mapping) {
        DEBUG("Processing read");
        for (const auto& path : path_extractor_(mapping)) {
            buffer_storages_[thread_index].AddPath(path, 1, false);
        }
        DEBUG("Read processed");
    }

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
    GappedPathExtractor(const Graph& g): g_(g), path_fixer_(g) {
    }

    vector<vector<EdgeId>> operator() (const MappingPath<EdgeId>& mapping) const {
        vector<EdgeId> corrected_path = path_fixer_.DeleteSameEdges(
                mapping.simple_path());
        corrected_path = FilterBadMappings(corrected_path, mapping);
        return FindReadPathWithGaps(mapping, corrected_path);
    }

private:

    size_t CountMappedEdgeSize(EdgeId edge, const MappingPath<EdgeId>& mapping_path, size_t& mapping_index) const {
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

    vector<EdgeId> FilterBadMappings(const vector<EdgeId>& corrected_path, const MappingPath<EdgeId>& mapping_path) const {
        vector<EdgeId> new_corrected_path;
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

    vector<vector<EdgeId>> FindReadPathWithGaps(const MappingPath<EdgeId>& mapping_path, vector<EdgeId>& corrected_path) const {
          if (mapping_path.size() == 0) {
              TRACE("read unmapped");
              return vector<vector<EdgeId>>();
          }
          vector<EdgeId> fixed_path = path_fixer_.TryFixPath(corrected_path);
          return SplitUnfixedPoints(fixed_path);
      }

    vector<vector<EdgeId>> SplitUnfixedPoints(vector<EdgeId>& path) const {
        vector<vector<EdgeId>> result;
        size_t prev_start = 0;
        for (size_t i = 1; i < path.size(); ++i) {
            if (g_.EdgeEnd(path[i - 1]) != g_.EdgeStart(path[i])) {
                    result.push_back(vector<EdgeId>(path.begin() + prev_start, path.begin() + i));
                    prev_start = i;
            }
        }
        result.push_back(vector<EdgeId>(path.begin() + prev_start, path.end()));
        return result;
    }
};

typedef std::function<vector<vector<EdgeId>> (const MappingPath<EdgeId>&)> PathExtractionF;

inline PathExtractionF ChooseProperReadPathExtractor(const Graph& g, io::LibraryType lib_type) {
    if (lib_type == io::LibraryType::PathExtendContigs || lib_type == io::LibraryType::TSLReads
        || lib_type == io::LibraryType::TrustedContigs || lib_type == io::LibraryType::UntrustedContigs) {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return GappedPathExtractor(g)(mapping);
        };
    } else {
        return [&] (const MappingPath<EdgeId>& mapping) {
            return vector<vector<EdgeId>>{ReadPathFinder<Graph>(g).FindReadPath(mapping)};
        };
    }
}

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
