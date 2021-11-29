//***************************************************************************
//* Copyright (c) 2015-2021 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alignment/sequence_mapper_notifier.hpp"
#include "alignment/sequence_mapper_fwd.hpp"

#include "assembly_graph/paths/mapping_path.hpp"
#include "assembly_graph/core/graph.hpp"
#include "io/reads/paired_read.hpp"
#include "io/reads/read_stream_vector.hpp"

#include "utils/perf/timetracer.hpp"
#include "pipeline/partask_mpi.hpp"

#include <string>
#include <vector>

namespace debruijn_graph {
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

    class MapLibFuncMPI : public MapLibBase {
    public:
        void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedRead>& streams) const override {
            MapLibMPI(listeners, mapper, streams);
        }
        void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleRead>& streams) const override {
            MapLibMPI(listeners, mapper, streams);
        }
        void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::SingleReadSeq>& streams) const override {
            MapLibMPI(listeners, mapper, streams);
        }
        void operator() (const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<io::PairedReadSeq>& streams) const override {
            MapLibMPI(listeners, mapper, streams);
        }

    private:
        template<class ReadType>
        void MapLibMPI(const std::vector<SequenceMapperListener*>& listeners, const SequenceMapper<Graph>& mapper, io::ReadStreamList<ReadType>& streams) const {
            SequenceMapperNotifierMPI notifier;
            for (auto listener: listeners) {
                notifier.Subscribe(listener);
            }
            notifier.ProcessLibrary(streams, mapper);
        }
    };

} // namespace debruijn_graph
