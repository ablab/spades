#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {
    class ScaffolderAnalysisStage : public spades::AssemblyStage {

     public:

        ScaffolderAnalysisStage() :
            AssemblyStage("Scaffold graph polishing", "scaffolder_analysis") {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) override;

     private:
        DECL_LOGGER("ScaffolderAnalysisStage")
    };

}