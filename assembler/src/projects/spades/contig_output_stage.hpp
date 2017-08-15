//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {


class ContigOutput : public spades::AssemblyStage {
private:
    bool output_paths_;
    std::string contigs_name_;
public:
    ContigOutput(bool output_paths = true, std::string contigs_name = cfg::get().co.contigs_name)
        : AssemblyStage("Contig Output", "contig_output"), output_paths_(output_paths), contigs_name_(contigs_name) { }

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *);
};

}
