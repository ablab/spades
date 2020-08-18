//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "simplification.hpp"
#include "position_storage.hpp"
#include "projects/unitig_coverage/profile_storage.hpp"
#include "io/graph/gfa_writer.hpp"
#include "toolchain/utils.hpp"

#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <string>
#include <numeric>
#include <sys/types.h>

using namespace debruijn_graph;

struct gcfg {
    gcfg() : k(0), RL(0),
             save_gfa(false), save_gp(false),
             use_cov_ratios(false),
             nthreads(omp_get_max_threads() / 2 + 1) {}

    unsigned k;
    unsigned RL;
    std::string bin_cov_str;
    std::string graph;
    std::string tmpdir;
    std::string edge_profile_fn;
    std::string stop_codons_fn;
    std::string deadends_fn;
    std::string outfile;
    bool save_gfa;
    bool save_gp;
    bool use_cov_ratios;
    unsigned nthreads;
};

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph. In GFA (ending with .gfa) or prefix to SPAdes graph pack"),
      cfg.outfile << value("output prefix"),
      option("--gfa").set(cfg.save_gfa) % "produce GFA output (default: true)",
      option("--spades-gp").set(cfg.save_gp) % "produce output graph pack in SPAdes internal format (default: false). "
                                                      "Recommended if bulges are removed to improve further read mapping. "
                                                      "In case GFA output is required with graph pack specify '--gfa'",
      option("--use-cov-ratios").set(cfg.use_cov_ratios) % "enable procedures based on unitig coverage ratios (default: false)",
      (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (required("--read-length") & integer("value", cfg.RL)) % "read length",
      (option("-c", "--coverage") & value("coverage", cfg.bin_cov_str)) % "estimated average (k+1-mer) bin coverage (default: 0.) "
                                                                          "or 'auto' (works only with '-d/--dead-ends' provided)",
      (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
      (option("-p", "--profile") & value("file", cfg.edge_profile_fn)) % "file with edge coverage profiles across multiple samples",
      (option("-s", "--stop-codons") & value("file", cfg.stop_codons_fn)) % "file stop codon positions",
      (option("-d", "--dead-ends") & value("file", cfg.deadends_fn)) % "while processing a subgraph -- file listing edges which are dead-ends in the original graph",
      (option("--tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use (default: <output prefix>.tmp)"
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
}

size_t DetermineSampleCnt(const std::string &profile_fn) {
    std::ifstream is(profile_fn);
    std::string line;
    std::getline(is, line);
    CHECK_FATAL_ERROR(is, "I/O problem while reading " << profile_fn);
    std::istringstream ss(line);
    std::string token;
    size_t i = 0;
    while (ss >> token) {
        ++i;
    }
    VERIFY(i > 0);
    return i - 1;
}

static std::set<std::string> ReadDeadendNames(const std::string &deadends_fn) {
    std::set<std::string> deadend_names;
    std::ifstream is(deadends_fn);
    std::string s;
    while (is >> s) {
        deadend_names.insert(s);
    }
    return deadend_names;
}

static void FillFlankingFromAverage(const Graph &g, FlankingCoverage<Graph> &flanking_cov) {
    for (EdgeId e : g.edges()) {
        auto raw_flank = g.coverage(e) * double(std::min(g.length(e), flanking_cov.averaging_range()));
        flanking_cov.SetRawCoverage(e, unsigned(math::round(raw_flank)));
    }
}

static bool IsDeadEnd(const Graph &g, VertexId v) {
    return g.IncomingEdgeCount(v) * g.OutgoingEdgeCount(v) == 0;
}

//TODO improve
//TODO think about self-conjugate sources/sinks and other connections between RC subgrahps
//TODO check consistency between total sink/source coverage estimates and think about weird cases
static double DetermineAvgCoverage(const Graph &g, const std::set<EdgeId> &/*undeadends*/) {
    double sum = 0.;
    for (auto it = g.ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        if (IsDeadEnd(g, g.EdgeStart(e))) {
            sum += g.coverage(e);
        }
        if (IsDeadEnd(g, g.EdgeEnd(e))) {
            sum += g.coverage(e);
        }
    }
    return sum / 2.;
}

//TODO set up reasonable flanking range
int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();
    START_BANNER("SPAdes-based standalone graph simplifier");

    try {
        unsigned nthreads = cfg.nthreads;
        unsigned k = cfg.k;
        std::string tmpdir = cfg.tmpdir.empty() ? cfg.outfile + ".tmp" : cfg.tmpdir;

        fs::make_dir(tmpdir);

        INFO("K-mer length set to " << k);

        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        debruijn_graph::GraphPack gp(k, tmpdir, 0);
        const auto &graph = gp.get<Graph>();

        INFO("Loading de Bruijn graph from " << cfg.graph);
        omnigraph::GraphElementFinder<Graph> element_finder(graph);
        gp.get_mutable<KmerMapper<Graph>>().Attach();

        io::EdgeLabelHelper<Graph> label_helper(element_finder,
                                                toolchain::LoadGraph(gp, cfg.graph));

        //Refilling flanking coverage to get same behavior while working with gfa graphs
        auto &flanking_cov = gp.get_mutable<FlankingCoverage<Graph>>();
        FillFlankingFromAverage(graph, flanking_cov);
        VERIFY(flanking_cov.IsAttached());

        //Loading and tracking edges that connect the subgraph to the rest of the graph
        std::set<EdgeId> undeadends;
        auto & edge_qual = gp.get_mutable<EdgeQuality<Graph>>();
        if (!cfg.deadends_fn.empty()) {
            auto deadend_names = ReadDeadendNames(cfg.deadends_fn);
            for (auto it = graph.ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
                EdgeId e = *it;
                if (IsDeadEnd(graph, graph.EdgeStart(e)) || IsDeadEnd(graph, graph.EdgeEnd(e))) {
                    if (!deadend_names.count(label_helper.label(e))) {
                        undeadends.insert(e);
                        undeadends.insert(graph.conjugate(e));
                    }
                }
            }
            //Using "quality" to track the 'undeadends'
            //TODO can use position storage for that instead of hacking edge_qual
            for (EdgeId e : undeadends) {
                VERIFY(undeadends.count(graph.conjugate(e)));
                edge_qual.AddQuality(e, 1000.);
            }
            edge_qual.Attach();
        }

        //Loading unitig profiles and keeping storage consistent
        typedef debruijn_graph::coverage_profiles::EdgeProfileStorage ProfileStorage;
        std::unique_ptr<ProfileStorage> profile_storage;
        if (!cfg.edge_profile_fn.empty()) {
            INFO("Loading edge profiles from " << cfg.edge_profile_fn);
            fs::CheckFileExistenceFATAL(cfg.edge_profile_fn);
            size_t sample_cnt = DetermineSampleCnt(cfg.edge_profile_fn);
            INFO("Sample count determined as " << sample_cnt);
            profile_storage = std::make_unique<ProfileStorage>(graph, sample_cnt);
            std::ifstream is(cfg.edge_profile_fn);
            profile_storage->Load(is, label_helper);
        }

        //Loading unitigs which were predicted to contain stop codons and keeping storage consistent
        using debruijn::simplification::PositionStorage;
        std::unique_ptr<PositionStorage> stop_codons_storage;
        if (!cfg.stop_codons_fn.empty()) {
            INFO("Loading stop codon positions from " << cfg.stop_codons_fn);
            fs::CheckFileExistenceFATAL(cfg.stop_codons_fn);
            stop_codons_storage = std::make_unique<PositionStorage>(graph);
            std::ifstream is(cfg.stop_codons_fn);
            stop_codons_storage->Load(is, label_helper);
        }

        //Setting the estimated mean coverage
        double bin_cov = 0.;
        if (cfg.bin_cov_str.empty()) {
            INFO("Mean coverage option was not specified. Using default value of " << bin_cov);
        } else {
            INFO("Estimated mean coverage was specified as " << cfg.bin_cov_str)
            if (cfg.bin_cov_str == "auto") {
                INFO("Trying to determine from coverage of sources and sinks");
                CHECK_FATAL_ERROR(!cfg.deadends_fn.empty(),
                           "Deadends option (-d/--dead-ends) was not specified. "
                           "Can only determine coverage while working with subgraphs!");
                bin_cov = DetermineAvgCoverage(graph, undeadends);
            } else {
                bin_cov = std::stod(cfg.bin_cov_str);
            }
        }

        //TODO parameterize
        double ec_bound = std::max(2.5, bin_cov / 50.);
        INFO("Erroneous connection coverage threshold set at " << ec_bound);

        //Setting up some simplification parameters
        debruijn::simplification::SimplifInfoContainer simplif_info;
        simplif_info.set_main_iteration(true);
        simplif_info.set_read_length(cfg.RL);
        simplif_info.set_chunk_cnt(200); //TODO magic constant
        simplif_info.set_detected_mean_coverage(bin_cov);
        simplif_info.set_detected_coverage_bound(ec_bound);

        //Simplification call
        debruijn::simplification::Simplify(gp, simplif_info, cfg.use_cov_ratios);

        INFO("Saving graph to " << cfg.outfile);

        fs::make_dirs(fs::parent_path(cfg.outfile));

        if (cfg.save_gp) {
            INFO("Saving graph pack in SPAdes binary format");
            io::binary::BasePackIO().Save(cfg.outfile, gp);
        }

        if (cfg.save_gfa || !cfg.save_gp) {
            INFO("Saving GFA");
            std::ofstream os(cfg.outfile + ".gfa");
            gfa::GFAWriter writer(graph, os);
            writer.WriteSegmentsAndLinks();
        }

        //Saving all the additional storages if were loaded
        if (profile_storage) {
            INFO("Saving profile storage");
            std::ofstream os(cfg.outfile + ".tsv");
            profile_storage->Save(os);
        }

        if (stop_codons_storage) {
            INFO("Saving codons");
            std::ofstream os(cfg.outfile + ".stops");
            stop_codons_storage->Save(os);
        }

        if (!cfg.deadends_fn.empty()) {
            INFO("Saving deadends");
            VERIFY_MSG(edge_qual.IsAttached(), "Edge quality got detached");
            std::ofstream deadends_os(cfg.outfile + ".deadends");
            for (auto it = graph.ConstEdgeBegin(/*canonical_only*/true); !it.IsEnd(); ++it) {
                EdgeId e = *it;
                if (IsDeadEnd(graph, graph.EdgeStart(e)) || IsDeadEnd(graph, graph.EdgeEnd(e))) {
                    if (edge_qual.IsPositiveQuality(e))
                        continue;
                    deadends_os << graph.int_id(e) << "\n";
                }
            }
        }

        INFO("Finished");

    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
