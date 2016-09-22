#pragma once

#include "tslr_pe.hpp"
#include "utils/logger/logger.hpp"
#include "barcode_mapper.hpp"
#include "common/pipeline/stage.hpp"

using namespace tslr_resolver;

namespace spades {
    class BarcodeMapConstructionStage : public AssemblyStage {
    private:
        size_t k_;

    public:

        BarcodeMapConstructionStage(size_t k) :
                AssemblyStage("Barcode map construction", "barcode_map_construction"),
                k_(k) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Barcode map construction started...");
            graph_pack.barcode_mapper = std::make_shared<HeadTailBarcodeMapper<SimpleBarcodeEntry>>
                    (graph_pack.g, cfg::get().ts_res.edge_tail_len);
            graph_pack.EnsureIndex();
            graph_pack.EnsureBasicMapping();
            graph_pack.barcode_mapper->FillMap(graph_pack.index, graph_pack.kmer_mapper);
            INFO("Barcode map construction finished.");
            INFO("Average barcode coverage: " +
                         std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));
        }
        DECL_LOGGER("BarcodeMapConstrusctionStage")
    };

} //spades
