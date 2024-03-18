
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class SeriesAnalysis : public spades::AssemblyStage {

public:
    SeriesAnalysis() : AssemblyStage("Series Analysis", "series_analysis") { }
    void load(graph_pack::GraphPack &, const std::filesystem::path &, const char *) override { }
    void save(const graph_pack::GraphPack &, const std::filesystem::path &, const char *) const override { }
    void run(graph_pack::GraphPack &gp, const char *) override;
};

}
