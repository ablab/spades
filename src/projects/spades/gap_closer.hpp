//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "alignment/sequence_mapper_notifier.hpp"
#include "pipeline/stage.hpp"
#include "io/reads/io_helper.hpp"

namespace debruijn_graph {

class GapClosingBase {
  protected:
    size_t num_readers = 0;
    virtual void ProcessLibrary(SequenceMapperListener* listener, const SequenceMapper<Graph>& mapper, io::BinaryPairedStreams& paired_streams) = 0;
  public:
    void execute(graph_pack::GraphPack &gp, const char *);
};

class GapClosing : public GapClosingBase, public spades::AssemblyStage {
  protected:
    void ProcessLibrary(SequenceMapperListener* listener, const SequenceMapper<Graph>& mapper, io::BinaryPairedStreams& paired_streams) override {
        SequenceMapperNotifier notifier;
        notifier.Subscribe(listener);
        notifier.ProcessLibrary(paired_streams, mapper);
    }

  public:
    GapClosing(const char* id)
            : AssemblyStage("Gap Closer", id) {}

    void run(graph_pack::GraphPack &gp, const char*) override;
};

}
