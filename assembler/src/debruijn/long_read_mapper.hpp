/*
 * long_read_mapper.hpp
 *
 *  Created on: Jun 17, 2013
 *      Author: andrey
 */

#ifndef LONG_READ_MAPPER_HPP_
#define LONG_READ_MAPPER_HPP_

#include "long_read_storage.hpp"
#include "sequence_mapper_notifier.hpp"

namespace debruijn_graph {

class SimpleLongReadMapper: public SequenceMapperListener {

public:

    SimpleLongReadMapper(conj_graph_pack& gp, PathStorage<conj_graph_pack::graph_t>& storage): gp_(gp), storage_(storage) {
        mapper_ = MapperInstance(gp_);
    }

    virtual ~SimpleLongReadMapper() {

    }

    virtual void StartProcessLibrary(size_t threads_count) {
        for (size_t i = 0; i < threads_count; ++i) {
            buffer_storages_.push_back(new PathStorage<conj_graph_pack::graph_t>(gp_.g));
        }
    }

    virtual void StopProcessLibrary() {
        for (size_t i = 0; i < buffer_storages_.size(); ++i) {
            MergeBuffer(i);
            delete buffer_storages_[i];
        }
        buffer_storages_.clear();
    }

    virtual void MergeBuffer(size_t thread_index) {
        DEBUG("Merge buffer " << thread_index << " with size " << buffer_storages_[thread_index]->size());
        storage_.AddStorage(*buffer_storages_[thread_index]);
        buffer_storages_[thread_index]->Clear();
        DEBUG("Now size " << storage_.size());
    }

    virtual void ProcessPairedRead(size_t ,
                                   const MappingPath<EdgeId>& ,
                                   const MappingPath<EdgeId>& ,
                                   size_t ) {
        //nothing to do
    }

    virtual void ProcessSingleRead(size_t thread_index,
                                   const MappingPath<EdgeId>& read) {
        vector<EdgeId> path = ProcessSingleRead(read);
        buffer_storages_[thread_index]->AddPath(path, 1, true);
    }

    PathStorage<conj_graph_pack::graph_t>& GetPaths() {
        return storage_;
    }

private:
    vector<EdgeId> ProcessSingleRead(const MappingPath<EdgeId>& path) const {
        return mapper_->FindReadPath(path);
    }

    conj_graph_pack& gp_;
    PathStorage<conj_graph_pack::graph_t>& storage_;
    std::shared_ptr<const NewExtendedSequenceMapper<conj_graph_pack::graph_t,
                    conj_graph_pack::index_t> > mapper_;
    vector<PathStorage<conj_graph_pack::graph_t>*> buffer_storages_;
};

}/*longreads*/

#endif /* LONG_READ_MAPPER_HPP_ */
