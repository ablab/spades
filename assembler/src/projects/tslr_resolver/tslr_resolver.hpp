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
        std::string path_to_reference_;
        //const Config &config_;

    public:

        TslrResolverStage(size_t k, const std::string& output_file, const std::string& path_to_reference) :
                AssemblyStage("TSLR repeat resolver", "tslr_repeat_resolver"),
                k_(k), output_file_(output_file), path_to_reference_(path_to_reference) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Resolver started...");
            graph_pack.edge_pos.Attach();
            debruijn_graph::FillPos(graph_pack, path_to_reference_, "");
            INFO(graph_pack.barcode_mapper->size());
            graph_pack.barcode_mapper->
                SerializeOverallDistribution(cfg::get().output_dir + "bardistr");
            INFO("Average barcode coverage before trimming: " + std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));
            graph_pack.barcode_mapper-> FilterByAbundance(3);

            INFO("Average barcode coverage after trimming: " + std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));
            graph_pack.barcode_mapper ->
                SerializeOverallDistribution(cfg::get().output_dir + "bardistr_after");
                
            tslr_resolver::LaunchBarcodePE (graph_pack);
            INFO("Resolver finished!");
        }
        DECL_LOGGER("TSLRResolverStage")
    };
} //spades
