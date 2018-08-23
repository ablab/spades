//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace spades {

class ReadConversion : public AssemblyStage {
public:
    ReadConversion()
            : AssemblyStage("Binary Read Conversion", "read_conversion") { }

    void run(debruijn_graph::conj_graph_pack &gp, const char *) override;
    void load(debruijn_graph::conj_graph_pack &, const std::string &load_from, const char *prefix = nullptr) override;
    void save(const debruijn_graph::conj_graph_pack &, const std::string &save_to, const char *prefix = nullptr) const override;
};

}
