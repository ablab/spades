/*
 * sequence_mapper_notifier.hpp
 *
 *  Created on: Sep 3, 2013
 *      Author: ira
 */

#ifndef SEQUENCE_MAPPER_NOTIFIER_HPP_
#define SEQUENCE_MAPPER_NOTIFIER_HPP_
#include <cstdlib>
#include <vector>
#include "sequence_mapper.hpp"
#include "short_read_mapper.hpp"
#include "io/paired_read.hpp"
#include "graph_pack.hpp"

namespace debruijn_graph {
//todo think if we still need all this
class SequenceMapperListener {
public:
    virtual void StartProcessLibrary(size_t threads_count) = 0;
    virtual void StopProcessLibrary() = 0;
    virtual void ProcessPairedRead(size_t thread_index, const MappingPath<EdgeId>& read1, const MappingPath<EdgeId>& read2, size_t dist) = 0;
    virtual void ProcessSingleRead(size_t thread_index, const MappingPath<EdgeId>& read) = 0;
    virtual void MergeBuffer(size_t thread_index) = 0;
    virtual ~SequenceMapperListener(){
    }
};

class SequenceMapperNotifier {
public:
    typedef SequenceMapper<conj_graph_pack::graph_t> SequenceMapperT;
//    typedef std::shared_ptr<
//            const NewExtendedSequenceMapper<conj_graph_pack::graph_t, conj_graph_pack::index_t> > Mapper;

    SequenceMapperNotifier(const conj_graph_pack& gp, bool send_true_distance = true)
            : gp_(gp), send_true_distance_(send_true_distance) {
    }

    void Subscribe(size_t lib_index, SequenceMapperListener* listener) {
        while ((int)lib_index >= (int)listeners_.size() - 1) {
            std::vector<SequenceMapperListener*> vect;
            listeners_.push_back(vect);
        }
        listeners_[lib_index].push_back(listener);
    }

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        size_t lib_index, const SequenceMapperT& mapper, size_t threads_count) {
        streams.reset();
        NotifyStartProcessLibrary(lib_index, threads_count);

        size_t counter = 0, n = 15;
        size_t fmem = get_free_memory();
#       pragma omp parallel for num_threads(threads_count) shared(counter)
        for (size_t ithread = 0; ithread < threads_count; ++ithread) {
            size_t size = 0;
            size_t limit = 200000;
            ReadType r;
            auto& stream = streams[ithread];
            stream.reset();
            bool end_of_stream = false;
            while (!end_of_stream) {
                end_of_stream = stream.eof();
                while (!end_of_stream && size < limit) {
                    stream >> r;
                    ++size;
                    NotifyProcessRead(r, mapper, lib_index, ithread);
                    end_of_stream = stream.eof();
                    // Stop filling buffer if the amount of available is smaller
                    // than half of free memory.
                    if (10 * get_free_memory() / 4 < fmem &&
                        size > 10000)
                        break;
                }
#               pragma omp critical
                {
                    counter += size;
                    if (counter >> n) {
                        INFO("Processed " << counter << " reads");
                        n += 1;
                    }
                    size = 0;
                    NotifyMergeBuffer(lib_index, ithread);
                }
            }
        }
        INFO("Processed " << counter << " reads");
        NotifyStopProcessLibrary(lib_index);
    }

private:
    template<class ReadType>
    void NotifyProcessRead(const ReadType& r, const SequenceMapperT& mapper, size_t ilib, size_t ithread) const;

    void NotifyStartProcessLibrary(size_t ilib, size_t thread_count) const {
        for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
                ++ilistener) {
            listeners_[ilib][ilistener]->StartProcessLibrary(thread_count);
        }
    }

    void NotifyStopProcessLibrary(size_t ilib) const {
        for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
                ++ilistener) {
            listeners_[ilib][ilistener]->StopProcessLibrary();
        }
    }

    void NotifyMergeBuffer(size_t ilib, size_t ithread) const {
        for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
                ++ilistener) {
            listeners_[ilib][ilistener]->MergeBuffer(ithread);
        }
    }
    const conj_graph_pack& gp_;
    bool send_true_distance_;
    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::PairedReadSeq& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const {

    const Sequence read1 = r.first().sequence();
    const Sequence read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper.MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper.MapSequence(read2);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        if (send_true_distance_) {
            listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.distance());
        }
        else {
            TRACE("Dist: " << r.second().size() << " - " << r.insert_size() << " = " << r.second().size() - r.insert_size());
            listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.second().size() - r.insert_size());
        }
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path1);
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path2);
    }
}
//TODO:delete it
template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::PairedRead& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const {
    const Sequence read1 = r.first().sequence();
    const Sequence read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper.MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper.MapSequence(read2);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        if (send_true_distance_) {
            listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.distance());
        }
        else {
            TRACE("Dist: " << r.second().size() << " - " << r.insert_size() << " = " << r.second().size() - r.insert_size());
            listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.second().size() - r.insert_size());
        }
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path1);
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path2);
    }
}

template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::SingleReadSeq& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const {
    Sequence read = r.sequence();
    MappingPath<EdgeId> path = mapper.MapSequence(read);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path);
    }
}

//TODO:delete it
template<>
inline void SequenceMapperNotifier::NotifyProcessRead(const io::SingleRead& r,
                                               const SequenceMapperT& mapper,
                                               size_t ilib,
                                               size_t ithread) const {
    Sequence read = r.sequence();
    MappingPath<EdgeId> path = mapper.MapSequence(read);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path);
    }
}

} /*debruijn_graph*/


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
