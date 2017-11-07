#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {
class ScaffoldGraphConstructionStage : public spades::AssemblyStage {

 public:

    ScaffoldGraphConstructionStage() :
        AssemblyStage("Scaffold graph construction", "scaffold_graph_construction") {
    }

    void run(debruijn_graph::conj_graph_pack &graph_pack, const char *);
    DECL_LOGGER("ReadCloudStatisticsStage")
};

} //debruijn_graph

