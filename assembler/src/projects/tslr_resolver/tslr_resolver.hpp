#pragma once

#include <modules/pipeline/stage.hpp>
#include <modules/stages/construction.hpp>
#include <tslr_pe.hpp>
#include <barcode_mapper.hpp>

namespace spades {
    class TslrResolverStage : public AssemblyStage {
        // public:
        //     typedef debruijn_graph::debruijn_config::tenx_resolver Config;
    private:
        size_t k_;
        std::string output_file_;
        //const Config &config_;

    public:

        TslrResolverStage(size_t k, const std::string& output_file) :
                AssemblyStage("TsrlResolver", "TSLR repeat resolver"),
                k_(k), output_file_(output_file) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Resolver started...");
            tslr_resolver::LaunchBarcodePE (graph_pack);
            INFO("Resolver finished!");
        }
        DECL_LOGGER("TSLRResolverStage")
    };
} //spades