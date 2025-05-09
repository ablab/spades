//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "mpi/alignment/sequence_mapper_notifier_mpi.hpp"
#include "mpi/pipeline/mpi_stage.hpp"
#include "io/reads/io_helper.hpp"
#include "projects/spades/gap_closer.hpp"

namespace debruijn_graph {
class GapClosingMPI : public GapClosingBase, public spades_mpi::MPIAssemblyStage {
protected:
    void ProcessLibrary(SequenceMapperListener *listener, const SequenceMapper<Graph> &mapper,
                        io::BinaryPairedStreams &paired_streams) override {
        SequenceMapperNotifierMPI notifier;
        notifier.Subscribe(listener);
        notifier.ProcessLibrary(paired_streams, mapper);
    }

  public:
    GapClosingMPI(const char* id)
            : MPIAssemblyStage("Gap Closer (parmap)", id) {
        num_readers = partask::overall_num_threads();
    }

    void run(graph_pack::GraphPack &gp, const char* s) override {
        execute(gp, s);
    }
};

}
