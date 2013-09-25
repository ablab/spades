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
#include "io/paired_read.hpp"

namespace debruijn_graph {

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
    typedef std::shared_ptr<
            const NewExtendedSequenceMapper<conj_graph_pack::graph_t, conj_graph_pack::index_t> > Mapper;
    SequenceMapperNotifier(const conj_graph_pack& gp)
            : gp_(gp) {
    }

    void Subscribe(size_t lib_index, SequenceMapperListener* listener) {
        while ((int)lib_index >= (int)listeners_.size() - 1) {
            std::vector<SequenceMapperListener*> vect;
            listeners_.push_back(vect);
        }
        listeners_[lib_index].push_back(listener);
    }

    template<class Read>
    void ProcessLibrary(io::ReadStreamVector<io::IReader<Read> >& streams,
            size_t lib_index, size_t threads_count) {
        streams.release();
        INFO("Processing reads (takes a while) " <<  threads_count);
        NotifyStartProcessLibrary(lib_index, threads_count);
        const Mapper mapper = MapperInstance(gp_);
        size_t counter = 0;
        #pragma omp parallel num_threads(threads_count)
        {
            #pragma omp for reduction(+: counter)
            for (size_t ithread = 0; ithread < threads_count; ++ithread) {
                INFO("thread " << ithread << " streams size " <<streams.size());
                size_t size = 0;
                size_t limit = 1000000;
                Read r;
                io::IReader<Read>& stream = streams[ithread];
                stream.reset();
                bool end_of_stream = false;
                while (!end_of_stream) {
                    end_of_stream = stream.eof();
                    while (!end_of_stream && size < limit) {
                        stream >> r;
                        ++size;
                        counter++;
                        NotifyProcessRead(mapper, lib_index, ithread, r);
                        end_of_stream = stream.eof();
                        if (size % 1000000 == 0) {
                            INFO("process " << counter << " reads");
                        }
                    }
                    #pragma omp critical
                    {
                        size = 0;
                        NotifyMergeBuffer(lib_index, ithread);
                    }
                }
            }
        }
        NotifyStopProcessLibrary(lib_index);
    }

private:
    template<class Read>
    void NotifyProcessRead(const Mapper mapper, size_t ilib, size_t ithread,
                           const Read& r) const;

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
    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};

template<>
void SequenceMapperNotifier::NotifyProcessRead(const Mapper mapper,
                                               size_t ilib,
                                               size_t ithread,
                                               const io::PairedReadSeq& r) const {
    const Sequence read1 = r.first().sequence();
    const Sequence read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper->MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper->MapSequence(read2);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.distance());
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path1);
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path2);
    }
}
//TODO:delete it
template<>
void SequenceMapperNotifier::NotifyProcessRead(const Mapper mapper,
                                               size_t ilib,
                                               size_t ithread,
                                               const io::PairedRead& r) const {
    const Sequence read1 = r.first().sequence();
    const Sequence read2 = r.second().sequence();
    MappingPath<EdgeId> path1 = mapper->MapSequence(read1);
    MappingPath<EdgeId> path2 = mapper->MapSequence(read2);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessPairedRead(ithread, path1, path2, r.distance());
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path1);
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path2);
    }
}

template<>
void SequenceMapperNotifier::NotifyProcessRead(const Mapper mapper,
                                               size_t ilib,
                                               size_t ithread,
                                               const io::SingleReadSeq& r) const {
    Sequence read = r.sequence();
    MappingPath<EdgeId> path = mapper->MapSequence(read);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path);
    }
}

//TODO:delete it
template<>
void SequenceMapperNotifier::NotifyProcessRead(const Mapper mapper,
                                               size_t ilib,
                                               size_t ithread,
                                               const io::SingleRead& r) const {
    Sequence read = r.sequence();
    MappingPath<EdgeId> path = mapper->MapSequence(read);
    for (size_t ilistener = 0; ilistener < listeners_[ilib].size();
            ++ilistener) {
        listeners_[ilib][ilistener]->ProcessSingleRead(ithread, path);
    }
}

} /*debruijn_graph*/


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
