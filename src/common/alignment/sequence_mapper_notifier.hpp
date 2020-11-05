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
        VERIFY_MSG(false, "Method Serialize is not implemented. Using default realization.");
    }

    virtual void Deserialize(std::istream&) {
        VERIFY_MSG(false, "Method Deserialize is not implemented. Using default realization.");
    }

    virtual void MergeFromStream(std::istream&) {
        VERIFY_MSG(false, "Method MergeFromStream is not implemented. Using default realization.");
    }

    virtual ~SequenceMapperListener() {}
};

inline void PyramidMergeMPI(SequenceMapperListener &listener) {
    size_t mpi_size = partask::world_size();
    size_t mpi_rank = partask::world_rank();
    const size_t deadbeef = 0xDEADBEEF;

    for (size_t step = 1; step < mpi_size; step *= 2) {
        if ((mpi_rank % (2*step) == 0) && (mpi_rank + step < mpi_size)) {
            partask::InputMPIStream is(mpi_rank + step);
            size_t sz;
            io::binary::BinRead(is, sz);
            VERIFY_MSG(sz == deadbeef, "Listener type: " << typeid(listener).name());
            listener.MergeFromStream(is);
            io::binary::BinRead(is, sz);
            VERIFY_MSG(sz == deadbeef, "Listener type: " << typeid(listener).name());
        } else if (mpi_rank % (2*step) == step) {
            partask::OutputMPIStream os(mpi_rank - step);
            io::binary::BinWrite(os, deadbeef);
            listener.Serialize(os);
            io::binary::BinWrite(os, deadbeef);
        }
    }
}

class SequenceMapperNotifier {
    static constexpr size_t BUFFER_SIZE = 200000;
public:
    typedef SequenceMapper<Graph> SequenceMapperT;

    typedef std::vector<SequenceMapperListener*> ListenersContainer;

    SequenceMapperNotifier(size_t lib_count = 1);

    void Subscribe(SequenceMapperListener* listener, size_t lib_index = 0);

    template<class ReadType>
    void ProcessLibraryMPI(io::ReadStreamList<ReadType>& streams,
                           size_t lib_index, const SequenceMapperT& mapper, size_t threads_count = 0) {
        INFO("ProcessLibraryMPI started");
        // Select streams
        std::vector<size_t> chunks;
        size_t mpi_size = partask::world_size();
        size_t mpi_rank = partask::world_rank();
        for (size_t i = 0; i < streams.size(); ++i) {
            if (i % mpi_size == mpi_rank) {
                chunks.push_back(i);
            }
        }
        INFO("Selected streams: " << chunks);
        auto local_streams = partask::create_empty_stream_list<ReadType>(chunks.size());
        partask::swap_streams(streams, local_streams, chunks);

        // Run ProcessLibrary
        INFO("Running ProcessLibrary");
        ProcessLibrary(local_streams, lib_index, mapper, threads_count);
        INFO("ProcessLibrary done");

        // Swap streams back
        partask::swap_streams(streams, local_streams, chunks);

        INFO("Merging results...");
        for (const auto& listener : listeners_[lib_index]) {
            INFO("Merging listener " << typeid(*listener).name());
            PyramidMergeMPI(*listener);
        }
        INFO("Listeners merged");

        const size_t deadbeef = 0xDEADBEEF;
        if (mpi_size > 1) {
            INFO("Syncing listeners...");
            if (mpi_rank == 0) {
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

  private:
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

    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};

} // namespace debruijn_graph


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
