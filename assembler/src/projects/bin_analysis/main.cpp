//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "bin_refinement.hpp"
#include "projects/subgraph_extractor/subgraph_extraction.hpp"

#include "io/graph/gfa_reader.hpp"
#include "io/id_mapper.hpp"
#include "io/binary/graph_pack.hpp"

#include "utils/logger/log_writers.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <unordered_map>
#include <string>
#include <numeric>
#include <sys/types.h>
#include <sys/stat.h>

static void create_console_logger() {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

static bool ends_with(const std::string &s, const std::string &p) {
    if (s.size() < p.size())
        return false;

    return (s.compare(s.size() - p.size(), p.size(), p) == 0);
}

using namespace debruijn_graph;

static void PrintGraphInfo(Graph &g) {
    size_t sz = 0;
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it)
        sz += 1;

    INFO("Graph loaded. Total vertices: " << g.size() << " Total edges: " << sz);
}

static std::unique_ptr<io::IdMapper<std::string>> LoadGraph(conj_graph_pack &gp,
                                                            const std::string &filename) {
    std::unique_ptr<io::IdMapper<std::string>> id_mapper;
    if (ends_with(filename, ".gfa")) {
        id_mapper = std::make_unique<io::IdMapper<std::string>>();
        gfa::GFAReader gfa(filename);
        INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
        gfa.to_graph(gp.g, &(*id_mapper));
    } else {
        io::binary::BasePackIO<Graph>().Load(filename, gp);
    }
    PrintGraphInfo(gp.g);
    return id_mapper;
}

static void DrawNeighbourComponents(conj_graph_pack &gp, const std::string &bin_contigs,
                             const EdgeQuality<Graph> &annotation,
                             const std::string &pics_path, size_t edge_len_bound) {
    gp.edge_pos.Attach();
    visualization::position_filler::FillPos(gp, bin_contigs, "bin_", /*with rc*/true);

    typedef bin_refinement::DFS<UnorientedNeighbourIteratorFactory<Graph>> ComponentCollectingDFS;
    std::set<EdgeId> considered;
    INFO("Considering long annotated edges");
    for (EdgeId e : annotation.PositiveQualEdges()) {
        if (gp.g.length(e) >= edge_len_bound && !considered.count(e)) {
            INFO("Collecting component from edge " << gp.g.str(e));
            ComponentCollectingDFS component_collector(gp.g, e, edge_len_bound, /*vertex limit*/50);
            bool vlimit_not_exceeded = component_collector.Run();
            utils::insert_all(considered, component_collector.border_sources());
            if (!vlimit_not_exceeded) {
                continue;
            }

            DEBUG("Writing to dot");
            bin_refinement::DrawComponent(gp.g,
                                          component_collector.reached(), pics_path + std::to_string(gp.g.int_id(e)) + ".dot",
                                          annotation.PositiveQualEdges(), bin_refinement::EdgeSet(), &gp.edge_pos);
            DEBUG("Written");
        }
    }
    gp.edge_pos.clear();
    gp.edge_pos.Detach();
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
    conj_graph_pack gp(K, tmpdir, 1, ref);

    INFO("Loading de Bruijn graph from " << graph_path);
    gp.kmer_mapper.Attach();

    auto id_mapper = LoadGraph(gp, graph_path);

    EdgeQuality<Graph> annotation(gp.g);
    gp.EnsureBasicMapping();

	//auto p = fs::append_path(out_prefix, "graph_pack");
    //io::binary::FullPackIO<Graph>().Save(p, gp);

    double base_cov = bin_refinement::AnnotateEdges(gp, annotation, bin_contigs);

    VERIFY(annotation.IsAttached());

    if (ref.size() > 0) {
        gp.EnsureDebugInfo();
    } else {
        VERIFY(!gp.edge_qual.IsAttached());
    }

    const std::string pics_path = out_prefix + "_neighbourhoods/";
    //const std::string pics_path = "";
    if (!pics_path.empty()) {
        fs::make_dirs(pics_path);
        DrawNeighbourComponents(gp, bin_contigs, annotation, pics_path, /*edge_length_bound*/2500);
    }

    INFO("Annotating extra edges");
    bin_refinement::AnnotateExtraEdges(gp, annotation, base_cov, out_prefix);

    auto gc = GraphComponent<Graph>::FromEdges(gp.g, annotation.PositiveQualEdges());
    gc = subgraph_extraction::ComponentExpander(gp.g).Expand(gc);
    io::EdgeNamingF<Graph> naming_f = id_mapper ? io::MapNamingF<Graph>(*id_mapper)
                                                : io::IdNamingF<Graph>();

    subgraph_extraction::WriteComponentWithDeadends(gc, out_prefix, naming_f);
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
      (option("-r", "--ref") & value("file", cfg.reference)) % "fasta file with reference sequence (for benchmarking purposes)",
      (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use",
      (option("-tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use"
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

    create_console_logger();
    START_BANNER("Bin Refinement");

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
