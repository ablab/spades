#pragma once

#include <barcode_map_construction.hpp>
#include <tslr_resolver.hpp>

namespace spades {

    void run_tslr_resolver(const std::string& path_to_tslr_dataset, const std::string& path_to_reference) {
        INFO("Starting from stage " << cfg::get().entry_point.c_str());

        debruijn_graph::conj_graph_pack conj_gp(cfg::get().K,
                                                cfg::get().tmp_dir,
                                                cfg::get().ds.reads.lib_count(),
                                                cfg::get().ds.reference_genome,
                                                cfg::get().flanking_range,
                                                cfg::get().pos.max_mapping_gap,
                                                cfg::get().pos.max_gap_diff);
        StageManager manager({cfg::get().developer_mode,
                              cfg::get().load_from,
                              cfg::get().output_saves});
        manager.add(new debruijn_graph::Construction())
                .add(new BarcodeMapConstructionStage(cfg::get().K, path_to_tslr_dataset))
                .add(new TslrResolverStage(cfg::get().K, cfg::get().output_dir + "resolver_output.fasta", path_to_reference));
        INFO("Output directory: " << cfg::get().output_dir);
        conj_gp.kmer_mapper.Attach();
        conj_gp.edge_pos.Attach();

        manager.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("TSLR resolver finished.");
    }
} //spades
