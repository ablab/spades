//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class AssemblyGraphOutput : public spades::AssemblyStage {
public:
    AssemblyGraphOutput()
        : AssemblyStage("Assembly Graph Output", "assembly_graph_output") { }

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *);
};


class ContigOutput : public spades::AssemblyStage {
public:
    ContigOutput()
        : AssemblyStage("Contig Output", "contig_output") { }

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *);
};

}