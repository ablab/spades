//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bin_refinement.hpp"
#include "toolchain/subgraph_utils.hpp"
#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"

#include "utils/segfault_handler.hpp"
#include "io/reads/file_reader.hpp"
#include "visualization/position_filler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <unordered_map>
#include <string>
#include <numeric>
#include <sys/types.h>

using namespace debruijn_graph;

static void DrawNeighbourComponents(GraphPack &gp, const std::string &bin_contigs,
                             const EdgeQuality<Graph> &annotation,
                             const std::string &pics_path, size_t edge_len_bound) {
    auto &edge_pos = gp.get_mutable<EdgesPositionHandler<Graph>>();
    const auto &graph = gp.get<Graph>();
    edge_pos.Attach();
    visualization::position_filler::FillPos(gp, bin_contigs, "bin_", /*with rc*/true);

    typedef bin_refinement::DFS<UnorientedNeighbourIteratorFactory<Graph>> ComponentCollectingDFS;
    std::set<EdgeId> considered;
    INFO("Considering long annotated edges");
    for (EdgeId e : annotation.PositiveQualEdges()) {
        if (graph.length(e) >= edge_len_bound && !considered.count(e)) {
            INFO("Collecting component from edge " << graph.str(e));
            ComponentCollectingDFS component_collector(graph, e, edge_len_bound, /*vertex limit*/50);
            bool vlimit_not_exceeded = component_collector.Run();
            utils::insert_all(considered, component_collector.border_sources());
            if (!vlimit_not_exceeded) {
                continue;
            }

            DEBUG("Writing to dot");
            bin_refinement::DrawComponent(graph,
                                          component_collector.reached(), pics_path + std::to_string(graph.int_id(e)) + ".dot",
                                          annotation.PositiveQualEdges(), bin_refinement::EdgeSet(), &edge_pos);
            DEBUG("Written");
        }
    }
    edge_pos.clear();
    edge_pos.Detach();
}

static void Run(size_t K, const std::string &graph_path,
                const std::string &bin_contigs, const std::string &out_prefix,
                const std::string &tmpdir, const std::string &reference_path = "") {
    using namespace bin_refinement;

    std::vector<std::string> ref;
    if (!reference_path.empty()) {
        fs::CheckFileExistenceFATAL(reference_path);
        io::FileReadStream genome_stream(reference_path);
        while (!genome_stream.eof()) {
            io::SingleRead r;
            genome_stream >> r;
            ref.push_back(r.GetSequenceString());
        }
    }

    //Fixme check if 1 is a good value
    GraphPack gp(K, tmpdir, 1, ref);
    const auto &graph = gp.get<Graph>();

    INFO("Loading de Bruijn graph from " << graph_path);
    omnigraph::GraphElementFinder<Graph> element_finder(graph);
    gp.get_mutable<KmerMapper<Graph>>().Attach();
    io::EdgeLabelHelper<Graph> label_helper(element_finder,
                                            toolchain::LoadGraph(gp, graph_path));

    gp.EnsureBasicMapping();

    EdgeQuality<Graph> annotation(graph);
    double base_cov = bin_refinement::AnnotateEdges(gp, annotation, bin_contigs);

    VERIFY(annotation.IsAttached());

    if (ref.size() > 0) {
        gp.EnsureDebugInfo();
    } else {
        VERIFY(!gp.get<EdgeQuality<Graph>>().IsAttached());
    }

    const std::string pics_path = out_prefix + "_neighbourhoods/";
    //const std::string pics_path = "";
    if (!pics_path.empty()) {
        fs::make_dirs(pics_path);
        DrawNeighbourComponents(gp, bin_contigs, annotation, pics_path, /*edge_length_bound*/2500);
    }

    INFO("Annotating extra edges");
    bin_refinement::AnnotateExtraEdges(gp, annotation, base_cov, out_prefix);

    auto gc = GraphComponent<Graph>::FromEdges(graph, annotation.PositiveQualEdges());
    gc = toolchain::ComponentExpander(graph).Expand(gc);

    toolchain::WriteComponentWithDeadends(gc, out_prefix, label_helper.edge_naming_f());
}

struct gcfg {
    gcfg()
        : k(55), tmpdir("tmp"), out_prefix(""),
          nthreads(omp_get_max_threads() / 2 + 1)
    {}

    unsigned k;
    std::string graph;
    std::string core_contigs;
    std::string reference;
    std::string tmpdir;
    std::string out_prefix;
    unsigned nthreads;
};

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph (in binary or GFA)"),
      cfg.core_contigs << value("path to contigs attributed to the bin"),
      cfg.out_prefix << value("output path prefix"),
      (option("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (option("-r", "--reference") & value("file", cfg.reference)) % "fasta file with reference sequence (for benchmarking purposes)",
      (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("--tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}

int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();
    START_BANNER("Experimental MAG Refinement with the assembly graphs");

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir;

        fs::make_dir(tmpdir);

        INFO("K-mer length set to " << k);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        Run(k, cfg.graph, cfg.core_contigs, cfg.out_prefix, tmpdir, cfg.reference);
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
