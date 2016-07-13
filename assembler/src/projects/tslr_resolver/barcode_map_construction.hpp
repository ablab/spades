#pragma once

#include <modules/pipeline/stage.hpp>
#include <modules/stages/construction.hpp>
#include <tslr_pe.hpp>
#include <barcode_mapper.hpp>

namespace spades {
    class BarcodeMapConstructionStage : public AssemblyStage {
        // public:
        //     typedef debruijn_graph::debruijn_config::tenx_resolver Config;
    private:
        size_t k_;
        std::string output_file_;
        const std::string tslr_dataset_;
        //const Config &config_;

    public:

        BarcodeMapConstructionStage(size_t k, const std::string& tslr_dataset) :
                AssemblyStage("Barcode Map Construction", "Barcode Map Construction"),
                k_(k), tslr_dataset_(tslr_dataset) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Barcode map construction started...");
            graph_pack.barcode_mapper.FillMap(tslr_dataset_, true);
            INFO("Barcode map construction finished.");
            INFO("Average barcode coverage: " + std::to_string(graph_pack.barcode_mapper.AverageBarcodeCoverage()));
        }
        DECL_LOGGER("BarcodeMapConstructionStage")
    };

} //spades
