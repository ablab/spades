//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SEQUENCE_MAPPER_NOTIFIER_HPP_
#define SEQUENCE_MAPPER_NOTIFIER_HPP_

#include "sequence_mapper.hpp"

#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/reads/paired_read.hpp"
#include "io/reads/read_stream_vector.hpp"

#include "utils/perf/timetracer.hpp"

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
    typedef SequenceMapper<Graph> SequenceMapperT;

    typedef std::vector<SequenceMapperListener*> ListenersContainer;

    SequenceMapperNotifier(const GraphPack &gp, size_t lib_count);

    void Subscribe(size_t lib_index, SequenceMapperListener* listener);

    template<class ReadType>
    void ProcessLibrary(io::ReadStreamList<ReadType>& streams,
                        size_t lib_index, const SequenceMapperT& mapper, size_t threads_count = 0) {
        std::string lib_str = std::to_string(lib_index);
        TIME_TRACE_SCOPE("SequenceMapperNotifier::ProcessLibrary", lib_str);
        if (threads_count == 0)
            threads_count = streams.size();

        streams.reset();
        NotifyStartProcessLibrary(lib_index, threads_count);
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

        for (size_t i = 0; i < threads_count; ++i)
            NotifyMergeBuffer(lib_index, i);

        INFO("Total " << counter << " reads processed");
        NotifyStopProcessLibrary(lib_index);
    }

private:
    template<class ReadType>
    void NotifyProcessRead(const ReadType& r, const SequenceMapperT& mapper, size_t ilib, size_t ithread) const;

    void NotifyStartProcessLibrary(size_t ilib, size_t thread_count) const;

    void NotifyStopProcessLibrary(size_t ilib) const;

    void NotifyMergeBuffer(size_t ilib, size_t ithread) const;

    const GraphPack& gp_;

    std::vector<std::vector<SequenceMapperListener*> > listeners_;  //first vector's size = count libs
};

} // namespace debruijn_graph


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
