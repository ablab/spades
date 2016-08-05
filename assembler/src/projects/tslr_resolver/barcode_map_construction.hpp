#pragma once

#include <modules/pipeline/stage.hpp>
#include <tslr_pe.hpp>
#include <barcode_mapper.hpp>

namespace spades {
    class BarcodeMapConstructionStage : public AssemblyStage {
        // public:
        //     typedef debruijn_graph::debruijn_config::tenx_resolver Config;
    private:
        size_t k_;
        const std::string tslr_dataset_;
        //const Config &config_;

    public:

        BarcodeMapConstructionStage(size_t k, const std::string& tslr_dataset) :
                AssemblyStage("Barcode map construction", "barcode_map_construction"),
                k_(k), tslr_dataset_(tslr_dataset) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Barcode map construction started...");
            graph_pack.barcode_mapper = std::make_shared<tslr_resolver::TrimmableBarcodeMapper>
                    (graph_pack.g, tslr_resolver::MapperType::Trimmable);
            bool debug_construction = cfg::get().ts_res.debug_construction;  //FIXME Remove that later
            graph_pack.barcode_mapper->FillMap(tslr_dataset_, graph_pack.index,
                                               graph_pack.kmer_mapper, debug_construction);
            INFO("Barcode map construction finished.");
            INFO("Average barcode coverage: " + std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));
        }
        DECL_LOGGER("BarcodeMapConstrusctionStage")
    };

} //spades
