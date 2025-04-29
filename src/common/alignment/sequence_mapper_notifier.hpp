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
#include "utils/stl_utils.hpp"

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

    virtual const std::string name() const {
        return utils::type_name(typeid(*this).name());
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
    std::vector<ListenersContainer> listeners_;  //first vector's size = count libs
};

class MapLibBase {
public:
    virtual void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedRead>& streams) const = 0;
    virtual void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleRead>& streams) const = 0;
    virtual void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleReadSeq>& streams) const = 0;
    virtual void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedReadSeq>& streams) const = 0;

    template<class Streams>
    void operator() (SequenceMapperListener* listener, const SequenceMapper<Graph>& mapper, Streams& streams) const {
        this->operator() (std::vector<SequenceMapperListener*>(1, listener), mapper, streams);
    }
};

class MapLibFunc : public MapLibBase {
public:
    void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedRead>& streams) const override {
        MapLib(listeners, mapper, streams);
    }
    void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleRead>& streams) const override {
        MapLib(listeners, mapper, streams);
    }
    void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleReadSeq>& streams) const override {
        MapLib(listeners, mapper, streams);
    }
    void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedReadSeq>& streams) const override {
        MapLib(listeners, mapper, streams);
    }

private:
    template<class ReadType>
    void MapLib(const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<ReadType>& streams) const {
        SequenceMapperNotifier notifier;
        for (auto listener: listeners) {
            notifier.Subscribe(listener);
        }
        notifier.ProcessLibrary(streams, mapper);
    }
};

} // namespace debruijn_graph


#endif /* SEQUENCE_MAPPER_NOTIFIER_HPP_ */
