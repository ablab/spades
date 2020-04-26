#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class SeriesAnalysis : public spades::AssemblyStage {

public:
    SeriesAnalysis() : AssemblyStage("Series Analysis", "series_analysis") { }

    void load(GraphPack &, const std::string &, const char *) { }

    void save(const GraphPack &, const std::string &, const char *) const { }

    void run(GraphPack &gp, const char *);
};

}
