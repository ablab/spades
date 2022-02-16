#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class SeriesAnalysis : public spades::AssemblyStage {

public:
    SeriesAnalysis() : AssemblyStage("Series Analysis", "series_analysis") { }

    void load(graph_pack::GraphPack &, const std::string &, const char *) override { }

    void save(const graph_pack::GraphPack &, const std::string &, const char *) const override { }

    void run(graph_pack::GraphPack &gp, const char *) override;
};

}
