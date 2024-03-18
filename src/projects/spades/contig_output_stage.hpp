//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {


class ContigOutput : public spades::AssemblyStage {
public:
    enum class Kind {
        BinaryContigs,
        EdgeSequences,
        GFAGraph,
        FASTGGraph,
        FinalContigs,
        PlasmidContigs,
        Scaffolds
    };
    typedef std::map<Kind, std::string> OutputList;

    ContigOutput(OutputList list)
            : AssemblyStage("Contig Output", "contig_output"),
              outputs_(std::move(list)) {}

    void save(const graph_pack::GraphPack &, const std::filesystem::path &, const char *) const override { }
    void run(graph_pack::GraphPack &gp, const char *) override;

private:
    OutputList outputs_;
};

}
