#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {

class SeriesAnalysis : public spades::AssemblyStage {

public:
    SeriesAnalysis() : AssemblyStage("Series Analysis", "series_analysis") { }

    void load(conj_graph_pack &, const std::string &, const char *) { }

    void save(const conj_graph_pack &, const std::string &, const char *) const { }

    void run(conj_graph_pack &gp, const char *);
};

}
