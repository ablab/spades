#include "barcode_index_construction.hpp"
#include <vector>
#include <common/barcode_index/barcode_info_extractor.hpp>
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/distribution_extractor_helper.hpp"
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
        INFO("Barcode index construction started...");
        const auto& dataset_info = cfg::get().ds;
        if (not HasReadClouds(dataset_info)) {
            INFO("Read cloud libraries have not been found. Skipping barcode index construction.")
            return;
        }
        vector<io::SequencingLibrary<debruijn_graph::config::LibraryData>> libs_10x;
        for (const auto& lib: dataset_info.reads) {
            if (lib.type() == io::LibraryType::Clouds10x) {
                libs_10x.push_back(lib);
            }
        }
        graph_pack.EnsureIndex();
        graph_pack.EnsureBasicMapping();
        FrameMapperBuilder mapper_builder(graph_pack.g,
                                          cfg::get().ts_res.edge_tail_len,
                                          cfg::get().ts_res.frame_size);

//        size_t nthreads = cfg::get().max_threads;
        mapper_builder.FillMap(libs_10x, graph_pack.index, graph_pack.kmer_mapper);
        graph_pack.barcode_mapper_ptr = mapper_builder.GetMapper();
        INFO("Barcode index construction finished.");
        FrameBarcodeIndexInfoExtractor extractor(graph_pack.barcode_mapper_ptr, graph_pack.g);
        INFO("Average barcode coverage: " + std::to_string(extractor.AverageBarcodeCoverage()));
        INFO("Estimating read cloud statistics using long edges");
        path_extend::cluster_model::ClusterStatisticsExtractorHelper cluster_extractor_helper(graph_pack,
                                                                                              cfg::get().max_threads);
        auto cluster_statistics_extractor = cluster_extractor_helper.GetStatisticsExtractor();
        auto distribution_pack = cluster_statistics_extractor.GetDistributionPack();
        graph_pack.read_cloud_distribution_pack = distribution_pack;
    }
}