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
    string contig_name_prefix_;

public:
    ContigOutput(bool output_paths = true, bool preliminary = false, const string& contig_name_prefix = "")
        : AssemblyStage("Contig Output", preliminary ? "preliminary_contig_output" : "contig_output"),
          output_paths_(output_paths), contig_name_prefix_(contig_name_prefix) { }

    void run(conj_graph_pack &gp, const char *);

};

}