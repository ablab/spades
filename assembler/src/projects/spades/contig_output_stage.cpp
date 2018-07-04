//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/path_extend/pe_resolver.hpp"
#include "contig_output_stage.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"

namespace debruijn_graph {

vector<path_extend::PathsWriterT> CreatePathsWriters(const std::string &fn_base,
                                                     path_extend::FastgPathWriter &fastg_writer) {
    using namespace path_extend;
    vector<PathsWriterT> writers;

    writers.push_back(ContigWriter::BasicFastaWriter(fn_base + ".fasta"));
    INFO("Outputting FastG paths to " << fn_base << ".paths");
    writers.push_back([=](const ScaffoldStorage& scaffold_storage) {
        fastg_writer.WritePaths(scaffold_storage, fn_base + ".paths");
    });
    return writers;
}

template<class Graph>
io::EdgeNamingF<Graph> PlasmidNamingF(io::EdgeNamingF<Graph> naming_f,
                                      const ConnectedComponentCounter &cc_counter) {
    return [=, &cc_counter](const Graph &g, EdgeId e) {
        return io::AddComponentId(naming_f(g, e), cc_counter.GetComponent(e));
    };
}

void ContigOutput::run(conj_graph_pack &gp, const char*) {
    using namespace path_extend;
    auto output_dir = cfg::get().output_dir;

    std::string gfa_fn = output_dir + "assembly_graph_with_scaffolds.gfa";
    INFO("Writing GFA to " << gfa_fn);

    std::ofstream os(gfa_fn);
    GFAPathWriter gfa_writer(gp.g, os,
                             cfg::get().pd ?
                             PlasmidNamingF<Graph>(io::IdNamingF<Graph>(), gp.components) :
                             io::IdNamingF<Graph>());
    gfa_writer.WriteSegmentsAndLinks();

    OutputEdgeSequences(gp.g, output_dir + "before_rr");

    INFO("Outputting FastG graph to " << output_dir << "assembly_graph.fastg");
    std::string fastg_fn = output_dir + "assembly_graph.fastg";

    FastgPathWriter fastg_writer(gp.g,
                                 fastg_fn,
                                 cfg::get().pd ?
                                 PlasmidNamingF<Graph>(io::BasicNamingF<Graph>(), gp.components) :
                                 io::BasicNamingF<Graph>());
    fastg_writer.WriteSegmentsAndLinks();

    if (output_paths_ && gp.contig_paths.size() != 0) {
        auto name_generator = MakeContigNameGenerator(cfg::get().mode, gp);
        ContigWriter writer(gp.g, name_generator);

        bool output_broken_scaffolds = cfg::get().pe_params.param_set.scaffolder_options.enabled &&
            cfg::get().use_scaffolder &&
            cfg::get().co.obs_mode != config::output_broken_scaffolds::none;

        if (output_broken_scaffolds) {
            int min_overlap = int(gp.g.k());
            if (cfg::get().co.obs_mode == config::output_broken_scaffolds::break_all) {
                min_overlap = int(gp.g.k());
            } else if (cfg::get().co.obs_mode == config::output_broken_scaffolds::break_gaps) {
                min_overlap = 0;
            } else {
                WARN("Unsupported contig output mode");
            }

            ScaffoldBreaker breaker(min_overlap);
            PathContainer broken_scaffolds;
            breaker.Break(gp.contig_paths, broken_scaffolds);

            //FIXME don't we want to use FinalizePaths here?
            GraphCoverageMap cover_map(gp.g, broken_scaffolds, true);
            Deduplicate(gp.g, broken_scaffolds, cover_map,
                    /*min_edge_len*/0,
                    /*max_path_diff*/0);
            broken_scaffolds.FilterEmptyPaths();
            broken_scaffolds.SortByLength();

            writer.OutputPaths(broken_scaffolds,
                               CreatePathsWriters(output_dir + contigs_name_,
                                                  fastg_writer));
        }

        auto writers = CreatePathsWriters(output_dir + cfg::get().co.scaffolds_name, fastg_writer);
        writers.push_back([&](const ScaffoldStorage &storage) {
            gfa_writer.WritePaths(storage);
        });
        writer.OutputPaths(gp.contig_paths, writers);
    } else {
        //FIXME weird logic
        OutputEdgeSequences(gp.g, output_dir + "simplified_contigs");
        OutputEdgeSequences(gp.g, output_dir + contigs_name_);
    }
}

}
