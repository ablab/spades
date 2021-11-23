//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SEQUENCE_MAPPER_NOTIFIER_HPP_
#define SEQUENCE_MAPPER_NOTIFIER_HPP_

#include "sequence_mapper_fwd.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/mapping_path.hpp"
#include "io/reads/paired_read.hpp"
#include "io/reads/read_stream_vector.hpp"
#include "utils/perf/timetracer.hpp"
#include "pipeline/partask_mpi.hpp"

#include <string>
#include <vector>

namespace debruijn_graph {
//todo think if we still need all this
class SequenceMapperListener {
public:
    virtual void StartProcessLibrary(size_t /* threads_count */) {}
    virtual void StopProcessLibrary() {}

    //TODO: think about read ierarchy
    virtual void ProcessPairedRead(size_t /* thread_index */, const io::PairedRead&  /* pr */,
                                   const omnigraph::MappingPath<EdgeId>& /* read1 */, const omnigraph::MappingPath<EdgeId>& /* read2 */) {}
    virtual void ProcessPairedRead(size_t /* thread_index */, const io::PairedReadSeq& /* pr */,
                                   const omnigraph::MappingPath<EdgeId>& /* read1 */, const omnigraph::MappingPath<EdgeId>& /* read2 */) {}
    virtual void ProcessSingleRead(size_t /* thread_index */, const io::SingleRead& /* r */, const omnigraph::MappingPath<EdgeId>& /* read */) {}
    virtual void ProcessSingleRead(size_t /* thread_index */, const io::SingleReadSeq& /* r */, const omnigraph::MappingPath<EdgeId>& /* read */) {}

    virtual void MergeBuffer(size_t /* thread_index */) {}

    virtual void Serialize(std::ostream&) const {
        VERIFY_MSG(false, "Serialize() is not implemented");
    }

    virtual void Deserialize(std::istream&) {
        VERIFY_MSG(false, "Deserialize() is not implemented");
    }

    virtual void MergeFromStream(std::istream&) {
        VERIFY_MSG(false, "MergeFromStream() is not implemented");
    }

    virtual const char* name() const {
        return typeid(*this).name();
    }

    virtual ~SequenceMapperListener() {}
};

class SequenceMapperNotifier {
    static constexpr size_t BUFFER_SIZE = 200000;
public:
    typedef SequenceMapper<Graph> SequenceMapperT;

    typedef std::vector<SequenceMapperListener*> ListenersContainer;

    SequenceMapperNotifier(size_t lib_count = 1);

    void Subscribe(SequenceMapperListener* listener, size_t lib_index = 0);

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        const SequenceMapperT& mapper, size_t threads_count = 0) {
        return ProcessLibrary(streams, 0, mapper, threads_count);
    }
    
    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        size_t lib_index, const SequenceMapperT& mapper, size_t threads_count = 0) {
        INFO("Starting sequence mapping");
        std::string lib_str = std::to_string(lib_index);
        TIME_TRACE_SCOPE("SequenceMapperNotifier::ProcessLibrary", lib_str);
        if (threads_count == 0)
            threads_count = streams.size();

        streams.reset();
        NotifyStartProcessLibrary(lib_index, streams.size());
        size_t counter = 0, n = 15;

        #pragma omp parallel for num_threads(threads_count) shared(counter)
        for (size_t i = 0; i < streams.size(); ++i) {
            size_t size = 0;
            ReadType r;
            auto& stream = streams[i];
            while (!stream.eof()) {
                if (size == BUFFER_SIZE) {
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

        for (size_t i = 0; i < streams.size(); ++i)
            NotifyMergeBuffer(lib_index, i);

        streams.close();
        INFO("Total " << counter << " reads processed");
        NotifyStopProcessLibrary(lib_index);
    }

private:
    template<class ReadType>
    void NotifyProcessRead(const ReadType& r, const SequenceMapperT& mapper, size_t ilib, size_t ithread) const;

    void NotifyStartProcessLibrary(size_t ilib, size_t thread_count) const;

    void NotifyStopProcessLibrary(size_t ilib) const;

    void NotifyMergeBuffer(size_t ilib, size_t ithread) const;

protected:
    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};


class SequenceMapperNotifierMPI : public SequenceMapperNotifier {
    void PyramidMergeMPI(SequenceMapperListener &listener);

public:
    using SequenceMapperNotifier::SequenceMapperNotifier;

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        size_t lib_index, const SequenceMapperT& mapper, size_t threads_count = 0) {
        INFO("ProcessLibraryMPI started");
        // Select streams
        std::vector<size_t> chunks = partask::chunks_rr(streams.size());
        INFO("Selected streams: " << chunks);

        partask::execute_on_subset(streams, chunks,
                                   [&](io::ReadStreamList<ReadType>& local_streams) {
                                       // Run ProcessLibrary
                                       INFO("Running ProcessLibrary");
                                       SequenceMapperNotifier::ProcessLibrary(local_streams, lib_index, mapper, threads_count);
                                       INFO("ProcessLibrary done");
                                   });

        INFO("Merging results...");
        for (const auto& listener : listeners_[lib_index]) {
            INFO("Merging listener " << listener->name());
            PyramidMergeMPI(*listener);
        }
        INFO("Listeners merged");

        if (partask::world_size() > 1) {
            const size_t deadbeef = 0xDEADBEEF;
            INFO("Syncing listeners...");
            if (partask::master()) {
                partask::OutputMPIStreamBcast os(0);
                for (const auto& listener : listeners_[lib_index]) {
                    io::binary::BinWrite(os, deadbeef);
                    listener->Serialize(os);
                    io::binary::BinWrite(os, deadbeef);
                }
            } else {
                partask::InputMPIStreamBcast is(0);
                for (const auto& listener : listeners_[lib_index]) {
                    size_t sz;
                    io::binary::BinRead(is, sz);
                    VERIFY(sz == deadbeef);
                    listener->Deserialize(is);
                    io::binary::BinRead(is, sz);
                    VERIFY(sz == deadbeef);
                }
            }
            INFO("Listeners synced");
        }
    }

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        const SequenceMapperT& mapper, size_t threads_count = 0) {
        return ProcessLibrary(streams, 0, mapper, threads_count);
    }
};

} // namespace debruijn_graph


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
