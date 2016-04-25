//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * long_read_mapper.hpp
 *
 *  Created on: Jun 17, 2013
 *      Author: andrey
 */

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "assembly_graph/graph_alignment/long_read_storage.hpp"
#include "assembly_graph/graph_alignment/sequence_mapper_notifier.hpp"

namespace debruijn_graph {

class AbstractLongReadMapper: public SequenceMapperListener {
public:
    AbstractLongReadMapper(conj_graph_pack& gp, PathStorage<conj_graph_pack::graph_t>& storage)
            : gp_(gp), storage_(storage), path_finder_(gp_.g) {
    }

    void StartProcessLibrary(size_t threads_count) override {
        for (size_t i = 0; i < threads_count; ++i)
            buffer_storages_.emplace_back(gp_.g);
    }

    void StopProcessLibrary() override {
        for (size_t i = 0; i < buffer_storages_.size(); ++i) {
            MergeBuffer(i);
        }
        buffer_storages_.clear();
    }

    void MergeBuffer(size_t thread_index) override {
        DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index].size());
        storage_.AddStorage(buffer_storages_[thread_index]);
        buffer_storages_[thread_index].Clear();
        DEBUG("Now size " << storage_.size());
    }

    void ProcessPairedRead(size_t ,
                           const io::PairedReadSeq&,
                           const MappingPath<EdgeId>& ,
                           const MappingPath<EdgeId>&) override {
        //nothing to do
    }

    void ProcessPairedRead(size_t ,
                           const io::PairedRead&,
                           const MappingPath<EdgeId>& ,
                           const MappingPath<EdgeId>&) override {
        //nothing to do
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

    PathStorage<conj_graph_pack::graph_t>& GetPaths() {
        return storage_;
    }

private:

    virtual void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& read) = 0;

protected:
    conj_graph_pack& gp_;
    PathStorage<conj_graph_pack::graph_t>& storage_;
    ReadPathFinder<conj_graph_pack::graph_t> path_finder_;
    std::vector<PathStorage<conj_graph_pack::graph_t> > buffer_storages_;

};

class SimpleLongReadMapper: public AbstractLongReadMapper {
public:
    SimpleLongReadMapper(conj_graph_pack& gp, PathStorage<conj_graph_pack::graph_t>& storage)
            : AbstractLongReadMapper(gp, storage) {
    }

private:

    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& read) override {
        vector<EdgeId> path = path_finder_.FindReadPath(read);
        buffer_storages_[thread_index].AddPath(path, 1, false);
    }
};

class FirstPELongReadMapper : public AbstractLongReadMapper {
private:
    typedef MappingPathFixer<Graph> GraphMappingPathFixer;
    const GraphMappingPathFixer path_fixer_;

public:
    FirstPELongReadMapper(conj_graph_pack& gp, PathStorage<conj_graph_pack::graph_t>& storage)
            : AbstractLongReadMapper(gp, storage), path_fixer_(gp.g) {
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
                mapping_index = mapping_path.size() - 1;
                break;
            }
        }
        size_t end_idx = mapping_index;
        size_t total_len = 0;
        for(size_t i = start_idx; i <= end_idx; ++i) {
            total_len += mapping_path[i].second.initial_range.size();
        }

        return total_len;
    }

    void FilterBadMappings(vector<EdgeId>& corrected_path, const MappingPath<EdgeId>& mapping_path) const {
        vector<EdgeId> new_corrected_path;
        size_t mapping_index = 0;
        for(auto edge : corrected_path) {
            size_t mapping_size = CountMappedEdgeSize(edge, mapping_path, mapping_index);
            size_t edge_len =  gp_.g.length(edge);
            //VERIFY(edge_len >= mapping_size);
            if((double) mapping_size / (double) edge_len > 0.3) {
                new_corrected_path.push_back(edge);
            }
        }
        std::swap(new_corrected_path, corrected_path);
    }


    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& read) override {
        vector<EdgeId> corrected_path = path_fixer_.DeleteSameEdges(
                read.simple_path());
        FilterBadMappings(corrected_path, read);
        vector<vector<EdgeId>> paths = path_finder_.FindReadPathWithGaps(read, corrected_path);
        for(auto path : paths) {
            buffer_storages_[thread_index].AddPath(path, 1, false);
        }
    }
};


}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
