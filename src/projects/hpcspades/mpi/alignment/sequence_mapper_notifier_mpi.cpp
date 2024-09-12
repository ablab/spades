//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "sequence_mapper_notifier_mpi.hpp"

namespace debruijn_graph {
void SequenceMapperNotifierMPI::PyramidMergeMPI(SequenceMapperListener &listener) {
    size_t mpi_size = partask::world_size();
    size_t mpi_rank = partask::world_rank();
    const size_t deadbeef = 0xDEADBEEF;

    for (size_t step = 1; step < mpi_size; step *= 2) {
        if ((mpi_rank % (2 * step) == 0) && (mpi_rank + step < mpi_size)) {
            partask::InputMPIStream is(mpi_rank + step);
            size_t sz;
            io::binary::BinRead(is, sz);
            VERIFY_MSG(sz == deadbeef, "Listener type: " << typeid(listener).name());
            listener.MergeFromStream(is);
            io::binary::BinRead(is, sz);
            VERIFY_MSG(sz == deadbeef, "Listener type: " << typeid(listener).name());
        } else if (mpi_rank % (2 * step) == step) {
            partask::OutputMPIStream os(mpi_rank - step);
            io::binary::BinWrite(os, deadbeef);
            listener.Serialize(os);
            io::binary::BinWrite(os, deadbeef);
        }
    }
}


void SequenceMapperNotifierMPI::SyncListeners(ListenersContainer &listeners) {
    if (partask::world_size() == 1)
        return;

    const size_t deadbeef = 0xDEADBEEF;
    INFO("Syncing listeners...");
    if (partask::master()) {
        TIME_TRACE_SCOPE("sync listeners", "master");
        partask::OutputMPIStreamBcast os(0);
        for (const auto &listener: listeners) {
            io::binary::BinWrite(os, deadbeef);
            listener->Serialize(os);
            io::binary::BinWrite(os, deadbeef);
        }
    } else {
        TIME_TRACE_SCOPE("sync listeners", "worker");
        partask::InputMPIStreamBcast is(0);
        for (const auto &listener: listeners) {
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

} // namespace debruijn_graph
