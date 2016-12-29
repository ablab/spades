#pragma once

#include "tslr_pe.hpp"
#include "utils/logger/logger.hpp"
#include "barcode_mapper.hpp"
#include "common/pipeline/stage.hpp"

using namespace tslr_resolver;

namespace spades {
    class BarcodeMapConstructionStage : public AssemblyStage {
    private:
        size_t k_;

    public:

        BarcodeMapConstructionStage(size_t k) :
                AssemblyStage("Barcode map construction", "barcode_map_construction"),
                k_(k) {
        }

        void TestTSLRBarcodeConstruction(debruijn_graph::conj_graph_pack &graph_pack) {
            INFO("Starting construction test...");
            vector <HeadTailMapperBuilder<EdgeEntry>> mapper_builders;
            for (int i = 0; i < 2; ++i) {
                HeadTailMapperBuilder<EdgeEntry> mapper_builder(graph_pack.g,
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
                INFO(mapper -> AverageBarcodeCoverage())
            }
        }


        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Barcode map construction started...");
            graph_pack.EnsureIndex();
            graph_pack.EnsureBasicMapping();
            INFO(cfg::get().ts_res.library_type)
            auto lib_type = GetLibType(cfg::get().ts_res.library_type);
            HeadTailMapperBuilder<EdgeEntry> mapper_builder(graph_pack.g,
                                                                     cfg::get().ts_res.edge_tail_len);

            size_t nthreads = cfg::get().max_threads;
            mapper_builder.FillMap(lib_type, graph_pack.index, graph_pack.kmer_mapper, nthreads);
            graph_pack.barcode_mapper = mapper_builder.GetMapper();
            INFO("Barcode map construction finished.");
            INFO("Average barcode coverage: " +
                         std::to_string(graph_pack.barcode_mapper->AverageBarcodeCoverage()));
        }
        DECL_LOGGER("BarcodeMapConstrusctionStage")
    };

} //spades
