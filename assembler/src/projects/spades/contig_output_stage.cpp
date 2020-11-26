//***************************************************************************
//* Copyright (c) 2015-2017 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contig_output_stage.hpp"

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/paths/bidirectional_path_io/bidirectional_path_output.hpp"
#include "modules/path_extend/pe_resolver.hpp"
#include "modules/path_extend/pe_utils.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "utils/filesystem/path_helper.hpp"

#include <unordered_set>

namespace {
using namespace debruijn_graph;

static constexpr double LARGE_FRACTION = 0.8;

bool CheckUsedPath(const path_extend::BidirectionalPath &path, std::unordered_set<EdgeId> &used_edges) {
    const Graph& g = path.g();
    size_t used_len = 0;
    size_t total_len = 0;
    size_t path_len = path.Size();
    for (size_t i = 0; i < path_len; i++) {
        size_t cur_len = g.length(path.At(i));
        total_len += cur_len;
        if (used_edges.count(path.At(i))) {
            used_len += cur_len;
        }
    }

    for (size_t i = 0; i < path_len; i++) {
        used_edges.insert(path.At(i));
        used_edges.insert(g.conjugate(path.At(i)));
    }
//FIXME: constant
    return (math::ge((double)used_len, (double)total_len * LARGE_FRACTION));
}

path_extend::PathContainer GetCircularScaffolds(const path_extend::PathContainer &sc_storage, std::unordered_set<EdgeId> &used_edges, size_t min_circular_length) {
    path_extend::PathContainer res;
    INFO("Banned " << used_edges.size() <<" used edges");
    for (const auto &entry : sc_storage) {
        const path_extend::BidirectionalPath &path = *entry.first;

        if (!path.IsCircular() ||
            CheckUsedPath(path, used_edges) ||
            path.Length() < min_circular_length)
            continue;

        res.Create(entry.first);
    }
    
    INFO("Got " << res.size() << " circular scaffolds");
    return res;
}

path_extend::PathContainer GetTipScaffolds(const path_extend::PathContainer &sc_storage, const std::unordered_set<VertexId> &forbidden_vertices, std::unordered_set<EdgeId> &used_edges, size_t min_linear_length) {
    path_extend::PathContainer res;
    for (const auto &entry : sc_storage) {
        const path_extend::BidirectionalPath &path = *entry.first;
        if (path.Length() < min_linear_length ||
            CheckUsedPath(path, used_edges) ||
            !forbidden_vertices.count(path.g().EdgeStart(path.Front())) ||
            !forbidden_vertices.count(path.g().EdgeEnd(path.Back())))
            continue;
        
        res.Create(entry.first);
        
    }

    INFO("Got " << res.size() << " linear scaffolds");
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
                                                          boost::optional<path_extend::FastgPathWriter> fastg_writer = boost::none,
                                                          boost::optional<path_extend::GFAPathWriter> gfa_writer = boost::none) {
    using namespace path_extend;
    std::vector<path_extend::PathsWriterT> res;
    res.push_back([&](const ScaffoldStorage& scaffold_storage) {
                      std::string fn = fn_base + ".fasta";
                      INFO("Outputting contigs to " << fn);
                      ContigWriter::WriteScaffolds(scaffold_storage, fn);
                  });

    if (fastg_writer) {
        res.push_back([&](const ScaffoldStorage& scaffold_storage) {
                          INFO("Outputting FastG paths to " << fn_base << ".paths");
                          fastg_writer->WritePaths(scaffold_storage, fn_base + ".paths");
                      });
    }

    if (gfa_writer) {
        res.push_back([&](const ScaffoldStorage &storage) {
                          INFO("Populating GFA with scaffold paths");
                          gfa_writer->WritePaths(storage);
                      });
    }

    return res;
}

}

namespace debruijn_graph {

void ContigOutput::run(GraphPack &gp, const char*) {
    using namespace path_extend;
    auto output_dir = cfg::get().output_dir;
    const auto &graph = gp.get<Graph>();

    if (outputs_.count(Kind::BinaryContigs)) {
        std::string contigs_output_dir = fs::append_path(output_dir, outputs_[Kind::BinaryContigs]);
        fs::make_dir(contigs_output_dir);
        io::ReadConverter::ConvertEdgeSequencesToBinary(graph, contigs_output_dir, cfg::get().max_threads);
    }

    if (outputs_.count(Kind::EdgeSequences)) {
        OutputEdgeSequences(graph, fs::append_path(output_dir, outputs_[Kind::EdgeSequences]));
    }

    std::unique_ptr<std::ostream> gfa_os;
    boost::optional<GFAPathWriter> gfa_writer;

    const auto &components = gp.get<ConnectedComponentCounter>();
    if (outputs_.count(Kind::GFAGraph)) {
        io::EdgeNamingF<Graph> naming_f =
                config::PipelineHelper::IsPlasmidPipeline(cfg::get().mode) && components.IsFilled()?
                PlasmidNamingF<Graph>(io::IdNamingF<Graph>(), components) :
                io::IdNamingF<Graph>();
        std::string gfa_fn = fs::append_path(output_dir, outputs_[Kind::GFAGraph] + ".gfa");

        gfa_os.reset(new std::ofstream(gfa_fn));
        gfa_writer.emplace(graph, *gfa_os, naming_f);
        INFO("Writing GFA graph to " << gfa_fn);
        gfa_writer->WriteSegmentsAndLinks();
    }

    boost::optional<FastgPathWriter> fastg_writer;
    if (outputs_.count(Kind::FASTGGraph)) {
        io::EdgeNamingF<Graph> naming_f =
                config::PipelineHelper::IsPlasmidPipeline(cfg::get().mode) && components.IsFilled()?
                PlasmidNamingF<Graph>(io::BasicNamingF<Graph>(), components) :
                io::BasicNamingF<Graph>();

        std::string fastg_fn = fs::append_path(output_dir, outputs_[Kind::FASTGGraph] + ".fastg");

        fastg_writer.emplace(graph, fastg_fn, naming_f);
        INFO("Outputting FastG graph to " << fastg_fn);
        fastg_writer->WriteSegmentsAndLinks();
    }

    const auto &contig_paths = gp.get<path_extend::PathContainer>("exSPAnder paths");
    bool output_contig_paths =
            (outputs_.count(Kind::Scaffolds) || outputs_.count(Kind::FinalContigs) || outputs_.count(Kind::PlasmidContigs)) &&
            contig_paths.size();

    if (output_contig_paths) {
        ContigWriter writer(graph, MakeContigNameGenerator(cfg::get().mode, gp));

        bool output_broken_scaffolds = cfg::get().pe_params.param_set.scaffolder_options.enabled &&
                                       cfg::get().use_scaffolder &&
                                       cfg::get().co.obs_mode != config::output_broken_scaffolds::none &&
                                       (outputs_.count(Kind::FinalContigs) || outputs_.count(Kind::PlasmidContigs));

        if (output_broken_scaffolds) {
            int min_overlap = int(gp.k());
            switch (cfg::get().co.obs_mode) {
            default:
                WARN("Unsupported contig output mode");
                break;
            case config::output_broken_scaffolds::break_all:
                min_overlap = int(gp.k());
                break;
            case config::output_broken_scaffolds::break_gaps:
                min_overlap = 0;
                break;
            }

            INFO("Breaking scaffolds");
            ScaffoldBreaker breaker(min_overlap);
            PathContainer broken_scaffolds;
            breaker.Break(contig_paths, broken_scaffolds);

            //FIXME don't we want to use FinalizePaths here?
            GraphCoverageMap cover_map(graph, broken_scaffolds, /* subscribe */ true);
            Deduplicate(graph, broken_scaffolds, cover_map,
                            /* min_edge_len */ 0, /* max_path_diff */ 0);
            broken_scaffolds.FilterEmptyPaths();
            broken_scaffolds.SortByLength();

            if (outputs_.count(Kind::FinalContigs))
                writer.OutputPaths(broken_scaffolds,
                                   CreatePathsWriters(fs::append_path(output_dir, outputs_[Kind::FinalContigs]),
                                                      fastg_writer));

            if (outputs_.count(Kind::PlasmidContigs)) {
                using UsedEdges = omnigraph::SmartContainer<std::unordered_set<EdgeId>, Graph>;
                if (!gp.count<UsedEdges>("used_edges"))
                    gp.add("used_edges", UsedEdges(graph));
                PathContainer circulars = GetCircularScaffolds(broken_scaffolds, gp.get_mutable<UsedEdges>("used_edges"), cfg::get().pd->min_circular_length);
                writer.OutputPaths(circulars,
                                   CreatePathsWriters(fs::append_path(output_dir, outputs_[Kind::PlasmidContigs] + ".circular"),
                                                          fastg_writer));
                if (cfg::get().pd->output_linear) {
                    if (gp.count<PathContainer>("Plasmid paths")) {
                        PathContainer &pcres = gp.get_mutable<PathContainer>("Plasmid paths");
                        writer.OutputPaths(pcres,
                                           CreatePathsWriters(
                                                   fs::append_path(output_dir, outputs_[Kind::PlasmidContigs] + ".linearrepeat"),
                                                   fastg_writer));
                    }

                    using ForbiddenVertices = omnigraph::SmartContainer<std::unordered_set<VertexId>, Graph>;
                    PathContainer linears;
                    if (gp.count<ForbiddenVertices>("forbidden_vertices"))
                        linears = GetTipScaffolds(broken_scaffolds, gp.get<ForbiddenVertices>("forbidden_vertices"),
                                                  gp.get_mutable<UsedEdges>("used_edges"), cfg::get().pd->min_linear_length);
                    writer.OutputPaths(linears,
                                       CreatePathsWriters(
                                               fs::append_path(output_dir, outputs_[Kind::PlasmidContigs] + ".linears"),
                                               fastg_writer));
                }
            }
        }

        if (outputs_.count(Kind::Scaffolds))
            writer.OutputPaths(contig_paths,
                               CreatePathsWriters(fs::append_path(output_dir, outputs_[Kind::Scaffolds]),
                                                  fastg_writer, gfa_writer));
    } else if (outputs_.count(Kind::FinalContigs)) {
        OutputEdgeSequences(graph, fs::append_path(output_dir, outputs_[Kind::FinalContigs]));
    }
}

}
