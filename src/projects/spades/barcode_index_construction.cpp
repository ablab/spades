#include "barcode_index_construction.hpp"

#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor_helper.hpp"

#include <vector>

namespace debruijn_graph {
    //todo remove from here!

    bool HasReadClouds(config::dataset ds) {
        bool has_read_clouds = false;
        for (const auto& lib: ds.reads) {
            if (lib.type() == io::LibraryType::Clouds10x) {
                has_read_clouds = true;
            }
        }
        return has_read_clouds;
    }

    void BarcodeMapConstructionStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        using path_extend::read_cloud::fragment_statistics::ClusterStatisticsExtractorHelper;

        INFO("Barcode index construction started...");
        const auto& dataset_info = cfg::get().ds;
        if (not HasReadClouds(dataset_info)) {
            INFO("Read cloud libraries have not been found. Skipping barcode index construction.")
            return;
        }
        size_t num_threads = cfg::get().max_threads;
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            auto &lib = cfg::get_writable().ds.reads[i];
            if (lib.type() == io::LibraryType::Clouds10x) {
                graph_pack.EnsureIndex();
                graph_pack.EnsureBasicMapping();
                std::vector <io::SingleStreamPtr> reads;
                for (const auto &read: lib.reads()) {
                    auto stream = io::EasyStream(read, false);
                    reads.push_back(stream);
                }
                FrameMapperBuilder<Graph> mapper_builder(graph_pack.barcode_mapper,
                                                         cfg::get().pe_params.read_cloud.edge_tail_len,
                                                         cfg::get().pe_params.read_cloud.frame_size);
                mapper_builder.FillMap(reads, graph_pack.index, graph_pack.kmer_mapper);
                INFO("Barcode index construction finished.");
                FrameBarcodeIndexInfoExtractor extractor(graph_pack.barcode_mapper, graph_pack.g);
                size_t length_threshold = cfg::get().pe_params.read_cloud.long_edge_length_lower_bound;
                INFO("Average barcode coverage: " + std::to_string(extractor.AverageBarcodeCoverage(length_threshold)));
                ClusterStatisticsExtractorHelper cluster_extractor_helper(graph_pack.g, graph_pack.barcode_mapper,
                                                                          cfg::get().pe_params.read_cloud, num_threads);
                auto cluster_statistics_extractor = cluster_extractor_helper.GetStatisticsExtractor();
                auto distribution_pack = cluster_statistics_extractor.GetDistributionPack();
                lib.data().read_cloud_info.fragment_length_distribution = distribution_pack.length_distribution_;
            }
        }
//        graph_pack.read_cloud_distribution_pack = distribution_pack;
    }
}