//
// Created by itolstoganov on 13.09.17.
//

#pragma once

#include "pipeline/stage.hpp"

namespace debruijn_graph {
    class ScaffolderAnalysisStage : public spades::AssemblyStage {

     public:

        ScaffolderAnalysisStage() :
            AssemblyStage("Scaffolder analysis", "scaffolder_analysis") {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) override;
        DECL_LOGGER("ScaffolderAnalysisStage")
    };
}