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
        const std::string tslr_dataset_;
        //const Config &config_;

    public:

        TslrResolverStage(size_t k, const std::string& output_file, const std::string& tslr_dataset) :
                AssemblyStage("TsrlResolver", "TSLR repeat resolver"),
                k_(k), output_file_(output_file), tslr_dataset_(tslr_dataset) {
        }

        void run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
            INFO("Barcode map construction started...");
            tslr_resolver::BarcodeMapper <io::SingleRead> bmapper (graph_pack, tslr_dataset_);
            INFO("Barcode map construction finished.");
        }
        DECL_LOGGER("TSLR Resolver Stage")
    };


    void run_tslr_resolver(const std::string path_to_tslr_dataset) {
        INFO("TSLR resolver started");
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
                .add(new TslrResolverStage(cfg::get().K, cfg::get().output_dir + "resolver_output.fasta",
                path_to_tslr_dataset));
        INFO("Output directory: " << cfg::get().output_dir);
        conj_gp.kmer_mapper.Attach();
        conj_gp.edge_pos.Attach();
        manager.run(conj_gp, cfg::get().entry_point.c_str());
        INFO("TSLR resolver finished.");
    }
} //spades