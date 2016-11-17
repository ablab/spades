#pragma once

#include <tslr_pe.hpp>
#include <barcode_mapper.hpp>
#include "common/pipeline/stage.hpp"
#include "extension_chooser_statisticts.hpp"

namespace spades {
    class TslrResolverStage : public AssemblyStage {
    private:
        size_t k_;
        std::string output_file_;
        std::string path_to_reference_;

        void LaunchStatisticsCounter(const debruijn_graph::conj_graph_pack& graph_pack) {
            std::string statistics_directory = cfg::get().output_dir + "/stats";
            auto stat_collector = BarcodeStatisticsCollector(graph_pack.g, graph_pack.barcode_mapper);
            size_t length_bound = 1000;
            size_t distance_bound = cfg::get().ts_res.distance_bound;
            stat_collector.FillEdgeStatistics(length_bound, distance_bound);
            stat_collector.SerializeEdgeStatistics(statistics_directory + "/edge_statistics");
            stat_collector.FillTrimmingStatistics();
            stat_collector.SerializeTrimmingStatistics(statistics_directory + "/trimming_statistics");
            stat_collector.FillReadStatistics(cfg::get().ts_res.edge_tail_len);
            stat_collector.SerializeReadStatistics(statistics_directory + "/read_statistics");
            stat_collector.FillBarcodeStatistics();
            stat_collector.SerializeBarcodeStatistics(statistics_directory + "/barcode_statistics");
        }

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
            INFO("Trimming threshold: " << trimming_threshold)
            graph_pack.barcode_mapper->FilterByAbundance(trimming_threshold);
            graph_pack.barcode_mapper->
                    SerializeOverallDistribution(cfg::get().output_dir + "barcode_distribution_after_trimming");
            INFO("Average barcode coverage after trimming: " +
                 std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));

            if (cfg::get().ts_res.library_type == "tslr")
                tslr_resolver::LaunchBarcodePE(graph_pack);
            else {
                LaunchStatisticsCounter(graph_pack);
            }

            INFO("Resolver finished!");
        }

        DECL_LOGGER("TSLRResolverStage")


    };
} //spades
