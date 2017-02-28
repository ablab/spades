//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_output_stage.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

namespace debruijn_graph {

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    auto output_dir = cfg::get().output_dir + contig_name_prefix_;

    OutputContigs(gp.g, output_dir + "before_rr", false);
    OutputContigsToFASTG(gp.g, output_dir + "assembly_graph", gp.components);

    if (output_paths_ && gp.contig_paths.size() != 0) {
        DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
        DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);

        auto name_generator = path_extend::MakeContigNameGenerator(cfg::get().mode, gp);
        path_extend::ContigWriter writer(gp.g, constructor, gp.components, name_generator);

        bool output_broken_scaffolds = cfg::get().pe_params.param_set.scaffolder_options.enabled &&
            cfg::get().use_scaffolder &&
            cfg::get().co.obs_mode != config::output_broken_scaffolds::none;

        if (output_broken_scaffolds) {
            int min_gap = 0;
            if (cfg::get().co.obs_mode == config::output_broken_scaffolds::break_all) {
                min_gap = 1;
            } else if (cfg::get().co.obs_mode == config::output_broken_scaffolds::break_gaps) {
                min_gap = int(gp.g.k());
            } else {
                WARN("Unsupported contig output mode");
            }

            path_extend::ScaffoldBreaker breaker(min_gap);
            path_extend::PathContainer broken_scaffolds;
            breaker.Break(gp.contig_paths, broken_scaffolds);
            writer.OutputPaths(broken_scaffolds, output_dir + cfg::get().co.contigs_name);
        }

        writer.OutputPaths(gp.contig_paths, output_dir + cfg::get().co.scaffolds_name);

        OutputContigsToGFA(gp.g, gp.contig_paths, output_dir + "assembly_graph");
    } else {
        OutputContigs(gp.g, output_dir + "simplified_contigs", cfg::get().use_unipaths);
        OutputContigs(gp.g, output_dir + cfg::get().co.contigs_name, false);
    }
}

}
