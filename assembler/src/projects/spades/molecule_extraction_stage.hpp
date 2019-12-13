#pragma once

#include "utils/logger/logger.hpp"
#include "common/pipeline/stage.hpp"


namespace debruijn_graph {
    class MoleculeExtractionStage : public spades::AssemblyStage {

    public:

        MoleculeExtractionStage() :
                AssemblyStage("Molecule extraction stage", "molecule_extraction") {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *);

    private:
        DECL_LOGGER("MoleculeExtractionStage")
    };

} //debruijn_graph
