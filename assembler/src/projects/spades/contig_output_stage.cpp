//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_output_stage.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

namespace debruijn_graph {

vector<path_extend::PathsWriterT> CreatePathsWriters(const std::string &fn_base,
                                                     path_extend::FastgWriter<Graph> &fastg_writer) {
    using namespace path_extend;
    vector<PathsWriterT> writers;

    writers.push_back(ContigWriter::BasicFastaWriter(fn_base + ".fasta"));
    INFO("Outputting FastG paths to " << fn_base << ".paths");
    writers.push_back([&](const ScaffoldStorage& scaffold_storage) {
        fastg_writer.WritePaths(scaffold_storage, fn_base + ".paths");
    });
    return writers;
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    auto output_dir = cfg::get().output_dir;

    std::string gfa_fn = output_dir + "assembly_graph_with_scaffolds.gfa";
    INFO("Writing GFA to " << gfa_fn);

    std::ofstream os(gfa_fn);
    path_extend::GFAWriter<Graph> gfa_writer(gp.g, os);
    gfa_writer.WriteSegmentsAndLinks();

    OutputEdgeSequences(gp.g, output_dir + "before_rr");

    INFO("Outputting FastG graph to " << output_dir << "assembly_graph.fastg");
    std::string fastg_fn = output_dir + "assembly_graph.fastg";
    path_extend::FastgWriter<Graph> fastg_writer(gp.g, gp.components);
    fastg_writer.WriteSegmentsAndLinks(fastg_fn);

    if (output_paths_ && gp.contig_paths.size() != 0) {
        auto name_generator = path_extend::MakeContigNameGenerator(cfg::get().mode, gp);
        path_extend::ContigWriter writer(gp.g, name_generator);

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
            writer.OutputPaths(broken_scaffolds,
                               CreatePathsWriters(output_dir + cfg::get().co.contigs_name,
                                                  fastg_writer));
        }

        auto writers = CreatePathsWriters(output_dir + cfg::get().co.scaffolds_name, fastg_writer);
        writers.push_back([&](const path_extend::ScaffoldStorage &storage) {
            gfa_writer.WritePaths(storage);
        });
        writer.OutputPaths(gp.contig_paths, writers);
    } else {
        //FIXME weird logic
        OutputEdgeSequences(gp.g, output_dir + "simplified_contigs");
        OutputEdgeSequences(gp.g, output_dir + cfg::get().co.contigs_name);
    }
}

}
