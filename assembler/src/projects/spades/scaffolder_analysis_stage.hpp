//
// Created by itolstoganov on 13.09.17.
//

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

        void GetGraphStorageReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &small_scaffold_graph,
                                          const path_extend::scaffold_graph::ScaffoldGraph &large_scaffold_graph,
                                          const debruijn_graph::conj_graph_pack& graph_pack) const;

        void PrintScaffoldGraphReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                             const debruijn_graph::conj_graph_pack& graph_pack,
                                             size_t length_threshold) const;

        DECL_LOGGER("ScaffolderAnalysisStage")
    };

}