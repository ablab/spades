//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SEQUENCE_MAPPER_NOTIFIER_HPP_
#define SEQUENCE_MAPPER_NOTIFIER_HPP_

#include "sequence_mapper.hpp"
#include "io/reads/paired_read.hpp"
#include "io/reads/read_stream_vector.hpp"
#include "pipeline/graph_pack.hpp"
#include "common/utils/memory_limit.hpp"

#include <vector>
#include <cstdlib>

namespace debruijn_graph {
//todo think if we still need all this
class SequenceMapperListener {
public:
    virtual void StartProcessLibrary(size_t /* threads_count */) {}
    virtual void StopProcessLibrary() {}

    //TODO: think about read ierarchy
    virtual void ProcessPairedRead(size_t /* thread_index */, const io::PairedRead&  /* pr */,
                                   const MappingPath<EdgeId>& /* read1 */, const MappingPath<EdgeId>& /* read2 */) {}
    virtual void ProcessPairedRead(size_t /* thread_index */, const io::PairedReadSeq& /* pr */,
                                   const MappingPath<EdgeId>& /* read1 */, const MappingPath<EdgeId>& /* read2 */) {}
    virtual void ProcessSingleRead(size_t /* thread_index */, const io::SingleRead& /* r */, const MappingPath<EdgeId>& /* read */) {}
    virtual void ProcessSingleRead(size_t /* thread_index */, const io::SingleReadSeq& /* r */, const MappingPath<EdgeId>& /* read */) {}

    virtual void MergeBuffer(size_t /* thread_index */) {}
    
    virtual ~SequenceMapperListener() {}
};

class SequenceMapperNotifier {
    static constexpr size_t BUFFER_SIZE = 200000;
public:
    typedef SequenceMapper<conj_graph_pack::graph_t> SequenceMapperT;

    typedef std::vector<SequenceMapperListener*> ListenersContainer;

    SequenceMapperNotifier(const conj_graph_pack& gp, size_t lib_count)
            : gp_(gp), listeners_(lib_count) { }

    void Subscribe(size_t lib_index, SequenceMapperListener* listener) {
        VERIFY(lib_index < listeners_.size());
        listeners_[lib_index].push_back(listener);
    }

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        size_t lib_index, const SequenceMapperT& mapper, size_t threads_count = 0) {
        if (threads_count == 0)
            threads_count = streams.size();

        streams.reset();
        NotifyStartProcessLibrary(lib_index, threads_count);
        size_t counter = 0, n = 15;
        size_t fmem = utils::get_free_memory();

        #pragma omp parallel for num_threads(threads_count) shared(counter)
        for (size_t i = 0; i < streams.size(); ++i) {
            size_t size = 0;
            ReadType r;
            auto& stream = streams[i];
            while (!stream.eof()) {
                if (size == BUFFER_SIZE || 
                    // Stop filling buffer if the amount of available is smaller
                    // than half of free memory.
                    (10 * utils::get_free_memory() / 4 < fmem && size > 10000)) {
                    #pragma omp critical
                    {
                        counter += size;
                        if (counter >> n) {
                            INFO("Processed " << counter << " reads");
                            n += 1;
                        }
                        size = 0;
                        NotifyMergeBuffer(lib_index, i);
                    }
                }
                stream >> r;
                ++size;
                NotifyProcessRead(r, mapper, lib_index, i);
            }
            #pragma omp atomic
            counter += size;
        }

        for (size_t i = 0; i < threads_count; ++i)
            NotifyMergeBuffer(lib_index, i);

        INFO("Total " << counter << " reads processed");
        NotifyStopProcessLibrary(lib_index);
    }

private:
    template<class ReadType>
    void NotifyProcessRead(const ReadType& r, const SequenceMapperT& mapper, size_t ilib, size_t ithread) const;

    void NotifyStartProcessLibrary(size_t ilib, size_t thread_count) const {
        for (const auto& listener : listeners_[ilib])
            listener->StartProcessLibrary(thread_count);
    }

    void NotifyStopProcessLibrary(size_t ilib) const {
        for (const auto& listener : listeners_[ilib])
            listener->StopProcessLibrary();
    }

    void NotifyMergeBuffer(size_t ilib, size_t ithread) const {
        for (const auto& listener : listeners_[ilib])
            listener->MergeBuffer(ithread);
    }
    const conj_graph_pack& gp_;

    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::PairedReadSeq& r,
                                                      const SequenceMapperT& mapper,
                                                      size_t ilib,
                                                      size_t ithread) const {

    const Sequence& read1 = r.first().sequence();
    const Sequence& read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper.MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper.MapSequence(read2);
    for (const auto& listener : listeners_[ilib]) {
        listener->ProcessPairedRead(ithread, r, path1, path2);
        listener->ProcessSingleRead(ithread, r.first(), path1);
        listener->ProcessSingleRead(ithread, r.second(), path2);
    }
}

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::PairedRead& r,
                                                      const SequenceMapperT& mapper,
                                                      size_t ilib,
                                                      size_t ithread) const {
    MappingPath<EdgeId> path1 = mapper.MapRead(r.first());
    MappingPath<EdgeId> path2 = mapper.MapRead(r.second());
    for (const auto& listener : listeners_[ilib]) {
        listener->ProcessPairedRead(ithread, r, path1, path2);
        listener->ProcessSingleRead(ithread, r.first(), path1);
        listener->ProcessSingleRead(ithread, r.second(), path2);
    }
}

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::SingleReadSeq& r,
                                                      const SequenceMapperT& mapper,
                                                      size_t ilib,
                                                      size_t ithread) const {
    const Sequence& read = r.sequence();
    MappingPath<EdgeId> path = mapper.MapSequence(read);
    for (const auto& listener : listeners_[ilib])
        listener->ProcessSingleRead(ithread, r, path);
}

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::SingleRead& r,
                                                      const SequenceMapperT& mapper,
                                                      size_t ilib,
                                                      size_t ithread) const {
    MappingPath<EdgeId> path = mapper.MapRead(r);
    for (const auto& listener : listeners_[ilib])
        listener->ProcessSingleRead(ithread, r, path);
}

} /*debruijn_graph*/


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
