//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/alignment/long_read_storage.hpp"
#include "assembly_graph/core/graph.hpp"
#include "pipeline/stage.hpp"

namespace debruijn_graph {

namespace gap_closing {
class GapStorage;
};

void PacbioAlignLibrary(const Graph& g,
                        const io::SequencingLibrary<config::LibraryData>& lib,
                        PathStorage<Graph>& path_storage,
                        gap_closing::GapStorage& gap_storage,
                        size_t thread_cnt, const config::pacbio_processor &pb);


class HybridLibrariesAligning : public spades::AssemblyStage {
public:
    HybridLibrariesAligning()
            : AssemblyStage("Hybrid Aligning", "hybrid_aligning") {
    }

    void run(GraphPack &gp, const char*) override;
    DECL_LOGGER("HybridAligning");
};

}

