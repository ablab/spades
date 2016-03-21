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

#include "data_structures/assembly_graph/graph_alignment/long_read_storage.hpp"
#include "data_structures/assembly_graph/graph_alignment/sequence_mapper_notifier.hpp"

namespace debruijn_graph {

class SimpleLongReadMapper: public SequenceMapperListener {
public:
    SimpleLongReadMapper(conj_graph_pack& gp, PathStorage<conj_graph_pack::graph_t>& storage)
            : gp_(gp), storage_(storage), path_finder_(gp_.g) {
        mapper_ = MapperInstance(gp_);
    }

    virtual ~SimpleLongReadMapper() {}

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

    void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& read) {
        vector<vector<EdgeId>> paths = path_finder_.FindReadPathWithGaps(read);
        for(auto path : paths) {
            buffer_storages_[thread_index].AddPath(path, 1, false);
        }
    }

    conj_graph_pack& gp_;
    PathStorage<conj_graph_pack::graph_t>& storage_;
    std::shared_ptr<const NewExtendedSequenceMapper<conj_graph_pack::graph_t,
                    conj_graph_pack::index_t> > mapper_;
    ReadPathFinder<conj_graph_pack::graph_t> path_finder_;
    std::vector<PathStorage<conj_graph_pack::graph_t> > buffer_storages_;
};

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
