//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_output_stage.hpp"

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "modules/path_extend/pe_resolver.hpp"
#include "io/dataset_support/read_converter.hpp"

#include <unordered_set>

namespace {
using namespace debruijn_graph;

bool CheckCircularPath(const path_extend::BidirectionalPath* path) {
    return (path->Size() > 0 && path->g().EdgeStart(path->Front()) == path->g().EdgeEnd(path->Back()));
}

bool CheckUsedPath(const path_extend::BidirectionalPath* path, std::unordered_set<EdgeId> &used_edges) {
    const Graph& g = path->g();
    size_t used_len = 0;
    size_t total_len = 0;
    size_t path_len = path->Size();
    for (size_t i = 0; i < path_len; i++) {
        size_t cur_len = g.length(path->At(i));
        total_len += cur_len;
        if (used_edges.count(path->At(i))) {
            used_len += cur_len;
        }
    }

    for (size_t i = 0; i < path_len; i++) {
        used_edges.insert(path->At(i));
        used_edges.insert(g.conjugate(path->At(i)));
    }
//FIXME: constant
    return (math::ge((double)used_len, (double)total_len * ContigOutput::LARGE_FRACTION));
}

path_extend::PathContainer GetCircularScaffolds(const path_extend::PathContainer &sc_storage, std::unordered_set<EdgeId> &used_edges) {
    path_extend::PathContainer res;
    INFO("banned " << used_edges.size() <<" edges");
    for (auto it = sc_storage.begin(); it != sc_storage.end(); it++) {
//FIXME: constant
        if (CheckCircularPath(it->first) && !CheckUsedPath(it->first, used_edges) && it->first->Length() >= 500) {
            path_extend::BidirectionalPath *p = new path_extend::BidirectionalPath(*it->first);
            path_extend::BidirectionalPath *cp = new path_extend::BidirectionalPath(p->Conjugate());
            res.AddPair(p, cp);
        }
    }
    INFO("got circular scaffs");
    return res;
}

path_extend::PathContainer GetTipScaffolds(const path_extend::PathContainer &sc_storage, const std::unordered_set<VertexId> &forbidden_vertices) {
    path_extend::PathContainer res;
    for (auto it = sc_storage.begin(); it != sc_storage.end(); it++) {
//FIXME: constant
        if (it->first->Length() > 0 &&
            forbidden_vertices.count(it->first->g().EdgeStart(it->first->Front())) &&
            forbidden_vertices.count(it->first->g().EdgeEnd(it->first->Back()))) {
            path_extend::BidirectionalPath *p = new path_extend::BidirectionalPath(*it->first);
            path_extend::BidirectionalPath *cp = new path_extend::BidirectionalPath(p->Conjugate());
            res.AddPair(p, cp);
        }
    }
    INFO("got suspicious linear tips scaffs");
    return res;

}

template<class Graph>
io::EdgeNamingF<Graph> PlasmidNamingF(io::EdgeNamingF<Graph> naming_f,
                                      const ConnectedComponentCounter &cc_counter) {
    return [=, &cc_counter](const Graph &g, EdgeId e) {
        return io::AddComponentId(naming_f(g, e), cc_counter.GetComponent(e));
    };
}

std::vector<path_extend::PathsWriterT> CreatePathsWriters(const std::string &fn_base,
                                                          path_extend::FastgPathWriter &fastg_writer) {
    using namespace path_extend;
    std::vector<PathsWriterT> writers;

    writers.push_back(ContigWriter::BasicFastaWriter(fn_base + ".fasta"));
    writers.push_back([=](const ScaffoldStorage& scaffold_storage) {
        INFO("Outputting FastG paths to " << fn_base << ".paths");
        fastg_writer.WritePaths(scaffold_storage, fn_base + ".paths");
    });
    return writers;
}

}

namespace debruijn_graph {
void ContigOutput::run(conj_graph_pack &gp, const char*) {
    using namespace path_extend;
    auto output_dir = cfg::get().output_dir;

    if (!final_iteration_) {
        std::string contigs_output_dir = output_dir + "simplified_contigs";
        fs::make_dir(contigs_output_dir);
        io::ReadConverter::ConvertEdgeSequencesToBinary(gp.g, contigs_output_dir, cfg::get().max_threads);
        return;
    }

    OutputEdgeSequences(gp.g, output_dir + "before_rr");

    io::EdgeNamingF<Graph> naming_f =
            config::PipelineHelper::IsPlasmidPipeline(cfg::get().mode) && gp.components.IsFilled()?
            PlasmidNamingF<Graph>(io::IdNamingF<Graph>(), gp.components) :
            io::IdNamingF<Graph>();

    std::string gfa_fn = output_dir + "assembly_graph_with_scaffolds.gfa";
    INFO("Writing GFA graph to " << gfa_fn);

    std::ofstream os(gfa_fn);
    GFAPathWriter gfa_writer(gp.g, os, naming_f);
    gfa_writer.WriteSegmentsAndLinks();

    std::string fastg_fn = output_dir + "assembly_graph.fastg";
    INFO("Outputting FastG graph to " << fastg_fn);

    FastgPathWriter fastg_writer(gp.g, fastg_fn, naming_f);
    fastg_writer.WriteSegmentsAndLinks();

    if (output_paths_ && gp.contig_paths.size() != 0) {
        ContigWriter writer(gp.g, MakeContigNameGenerator(cfg::get().mode, gp));

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
            if (cfg::get().mode == config::pipeline_type::metaplasmid) {
                using UsedEdges = SmartContainer<std::unordered_set<EdgeId>, Graph>;
                if (!gp.count<UsedEdges>("used_edges")) {
                    gp.add("used_edges", UsedEdges(gp.g));
                }
                PathContainer circulars = GetCircularScaffolds(broken_scaffolds, gp.get_mutable<UsedEdges>("used_edges"));
                writer.OutputPaths(circulars,
                                   CreatePathsWriters(output_dir + contigs_name_ + ".circular",
                                                      fastg_writer));
                using ForbiddenVertices = SmartContainer<std::unordered_set<VertexId>, Graph>;
                PathContainer linears;
                if (gp.count<ForbiddenVertices>("forbidden_vertices"))
                    linears = GetTipScaffolds(broken_scaffolds, gp.get_const<ForbiddenVertices>("forbidden_vertices"));
                writer.OutputPaths(linears,
                                   CreatePathsWriters(output_dir + contigs_name_ + ".linears",
                                                      fastg_writer));
            }
        }

        auto writers = CreatePathsWriters(output_dir + cfg::get().co.scaffolds_name, fastg_writer);
        writers.push_back([&](const ScaffoldStorage &storage) {
            gfa_writer.WritePaths(storage);
        });
        writer.OutputPaths(gp.contig_paths, writers);

    } else {
        OutputEdgeSequences(gp.g, output_dir + contigs_name_);
    }
}

}
