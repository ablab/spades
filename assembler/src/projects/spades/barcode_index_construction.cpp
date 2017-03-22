#include "barcode_index_construction.hpp"
#include <vector>
#include <common/barcode_index/barcode_info_extractor.hpp>

namespace debruijn_graph {
    //todo remove from here!
    void TestTSLRBarcodeConstruction(debruijn_graph::conj_graph_pack &graph_pack) {
        INFO("Starting construction test...");
        vector <SimpleMapperBuilder> mapper_builders;
        for (int i = 0; i < 2; ++i) {
            SimpleMapperBuilder mapper_builder(graph_pack.g,
                                                            cfg::get().ts_res.edge_tail_len);
            mapper_builders.push_back(mapper_builder);
        }
        mapper_builders[0].FillMapFromDemultiplexedDataset(graph_pack.index, graph_pack.kmer_mapper);
        mapper_builders[1].FillMapUsingSubIndex(graph_pack.index, graph_pack.kmer_mapper);
        mapper_builders[2].FillMapUsingKmerMultisetParallel(graph_pack.index,
                                                            graph_pack.kmer_mapper,
                                                            cfg::get().max_threads);
        for (int i = 0; i < 2; ++i) {
            auto mapper = mapper_builders[i].GetMapper();
            INFO(std::to_string(i) + ":");
        }
    }

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
        vector<io::SequencingLibrary<debruijn_graph::config::DataSetData>> libs_10x;
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
    }
}