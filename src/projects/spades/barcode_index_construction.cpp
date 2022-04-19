#include "barcode_index_construction.hpp"

#include "barcode_index/barcode_info_extractor.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor_helper.hpp"

#include <vector>

namespace debruijn_graph {
    //todo remove from here!

    size_t ReadCloudLibsCount(config::dataset ds) {
        size_t result = 0;
        for (const auto& lib: ds.reads) {
            if (lib.type() == io::LibraryType::Clouds10x) {
                ++result;
            }
        }
        return result;
    }

    void BarcodeMapConstructionStage::run(debruijn_graph::GraphPack &gp, const char *) {
        using path_extend::read_cloud::fragment_statistics::ClusterStatisticsExtractorHelper;

        INFO("Barcode index construction started...");
        const auto& dataset_info = cfg::get().ds;
        size_t cloud_libs_count = ReadCloudLibsCount(dataset_info);
        if (cloud_libs_count == 0) {
            INFO("Read cloud libraries have not been found. Skipping barcode index construction.")
            return;
        } else if (cloud_libs_count > 1) {
            INFO("Multiple read cloud libraries are currently not supported. Skipping barcode index construction");
            return;
        }
        size_t num_threads = cfg::get().max_threads;
        size_t frame_size = cfg::get().pe_params.read_cloud.frame_size;
        const std::vector<string> barcode_prefices = {"BC:Z:", "BX:Z:"};

        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            auto &lib = cfg::get_writable().ds.reads[i];
            if (lib.type() == io::LibraryType::Clouds10x) {
                gp.EnsureBasicMapping();
                FrameConcurrentBarcodeIndexBuffer<Graph> buffer(gp.get<Graph>(), frame_size);
                ConcurrentBufferFiller buffer_filler(gp.get<Graph>(), buffer, *MapperInstance(gp), barcode_prefices, 0);
                FrameBarcodeIndexBuilder barcode_index_builder(gp.get<Graph>(), *MapperInstance(gp), barcode_prefices, frame_size, num_threads);
                auto &barcode_mapper = gp.get_mutable<barcode_index::FrameBarcodeIndex<Graph>>();
                barcode_index_builder.ConstructBarcodeIndex(barcode_mapper, lib);
                INFO("Barcode index construction finished.");
                FrameBarcodeIndexInfoExtractor extractor(barcode_mapper, gp.get<Graph>());
                size_t length_threshold = cfg::get().pe_params.read_cloud.long_edge_length_lower_bound;
                INFO("Average barcode coverage: " + std::to_string(extractor.AverageBarcodeCoverage(length_threshold)));
                ClusterStatisticsExtractorHelper cluster_extractor_helper(gp.get<Graph>(), barcode_mapper,
                                                                          cfg::get().pe_params.read_cloud, num_threads);
                auto cluster_statistics_extractor = cluster_extractor_helper.GetStatisticsExtractor();
                auto distribution_pack = cluster_statistics_extractor.GetDistributionPack();
                lib.data().read_cloud_info.fragment_length_distribution = distribution_pack.length_distribution_;
            }
        }
//        graph_pack.read_cloud_distribution_pack = distribution_pack;
    }
}