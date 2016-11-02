#pragma once

#include <tslr_pe.hpp>
#include <barcode_mapper.hpp>
#include "common/pipeline/stage.hpp"

namespace spades {
    class TslrResolverStage : public AssemblyStage {
    private:
        size_t k_;
        std::string output_file_;
        std::string path_to_reference_;

    public:

        TslrResolverStage(size_t k, const std::string &output_file) :
                AssemblyStage("TSLR repeat resolver", "tslr_repeat_resolver"),
                k_(k), output_file_(output_file) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Resolver started...");
            graph_pack.barcode_mapper->
                    SerializeOverallDistribution(cfg::get().output_dir + "barcode_distribution");
            INFO("Average barcode coverage before trimming: " +
                 std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));

            //todo get trimming threshold from distribution
            size_t trimming_threshold = cfg::get().ts_res.trimming_threshold;
            graph_pack.barcode_mapper->FilterByAbundance(trimming_threshold);


            tslr_resolver::LaunchBarcodePE(graph_pack);
            INFO("Resolver finished!");
        }

        DECL_LOGGER("TSLRResolverStage")


    };
} //spades
