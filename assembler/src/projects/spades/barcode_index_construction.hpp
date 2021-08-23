#pragma once

#include "common/barcode_index/barcode_index_builder.hpp"
#include "utils/logger/logger.hpp"
#include "common/pipeline/stage.hpp"

using namespace barcode_index;

namespace debruijn_graph {
    class BarcodeMapConstructionStage : public spades::AssemblyStage {

    public:

        BarcodeMapConstructionStage() :
                AssemblyStage("Barcode map construction", "barcode_map_construction") {
        }

        void run(debruijn_graph::GraphPack &gp, const char *);
        DECL_LOGGER("BarcodeMapConstrusctionStage")
    };

} //debruijn_graph
