//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "projects/edge_profiles/profile_storage.hpp"
#include "io/dataset_support/dataset_readers.hpp"
#include "io/graph/gfa_writer.hpp"
#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"
#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "pipeline/config_struct.hpp"

#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <string>
#include <numeric>
#include <sys/types.h>
#include <sys/stat.h>

using namespace debruijn_graph;

typedef config::debruijn_config::simplification::bulge_remover BRConfig;
typedef config::debruijn_config::simplification::relative_coverage_comp_remover RCCConfig;
typedef config::debruijn_config::simplification::relative_coverage_edge_disconnector REDConfig;

static BRConfig default_br_config(bool enabled = false) {
    BRConfig config;
    config.enabled = enabled;
    config.main_iteration_only = false;
    config.max_bulge_length_coefficient = 4;
    config.max_additive_length_coefficient = 0;
    config.max_coverage = 1000.;
    config.max_relative_coverage = 1.2;
    config.max_delta = 3;
    config.max_number_edges = std::numeric_limits<size_t>::max();
    config.dijkstra_vertex_limit = std::numeric_limits<size_t>::max();
    config.max_relative_delta = 0.1;
    config.parallel = false;
    config.buff_size = 10000;
    config.buff_cov_diff = 2.;
    config.buff_cov_rel_diff = 0.2;
    return config;
}

static RCCConfig default_rcc_config(bool enabled = false) {
    RCCConfig rcc_config;
    rcc_config.enabled = enabled;
    rcc_config.coverage_gap = 50.;
    rcc_config.length_coeff = 3.0;
    rcc_config.tip_allowing_length_coeff = 5.0;
    rcc_config.vertex_count_limit = 100;
    rcc_config.max_ec_length_coefficient = 300;
    rcc_config.max_coverage_coeff = -1.0;
    return rcc_config;
}

static REDConfig default_red_config(bool enabled = false) {
    REDConfig red_config;
    red_config.enabled = enabled;
    red_config.diff_mult = 75.;
    red_config.edge_sum = 10000;
    red_config.unconditional_diff_mult = 100.;
    return red_config;
}

//enum class output_type {
//    unitigs, fastg, gfa, spades, spades_pack
//};

struct gcfg {
    gcfg()
        : k(0), RL(0), save_gfa(false), rel_cov_proc_enabled(false),
          nthreads(omp_get_max_threads() / 2 + 1)
    {}

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
    bool rel_cov_proc_enabled;
    unsigned nthreads;
//    output_type mode;
};

//TODO normalize -/-- usage
static void process_cmdline(int argc, char **argv, gcfg &cfg) {
  using namespace clipp;

  auto cli = (
      cfg.graph << value("graph. In GFA (ending with .gfa) or prefix to SPAdes graph pack"),
      //cfg.outfile << value("output filename/prefix (in case of --spades-gp)"),
      cfg.outfile << value("output prefix"),
      option("--gfa").set(cfg.save_gfa) % "produce GFA output (default: false)",
      option("--rel-cov-proc").set(cfg.rel_cov_proc_enabled) % "enable procedures based on unitig coverage ratios (default: false)",
      (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
      (required("-read-length") & integer("value", cfg.RL)) % "read length",
      (option("-c", "--coverage") & value("coverage", cfg.bin_cov_str)) % "estimated average (k+1-mer) bin coverage (default: 0.) "
                                                                          "or 'auto' (works only with '-d/--dead-ends' provided)",
      (option("-t") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
      (option("-p", "--profile") & value("file", cfg.edge_profile_fn)) % "file with edge coverage profiles across multiple samples",
      (option("-s", "--stop-codons") & value("file", cfg.stop_codons_fn)) % "file stop codon positions",
      (option("-d", "--dead-ends") & value("file", cfg.deadends_fn)) % "while processing a subgraph -- file listing edges which are dead-ends in the original graph",
      (option("-tmpdir") & value("dir", cfg.tmpdir)) % "scratch directory to use (default: <output prefix>.tmp)"
//      one_of(option("--unitigs").set(cfg.mode, output_type::unitigs) % "produce unitigs (default)",
//             option("--fastg").set(cfg.mode, output_type::fastg) % "produce graph in FASTG format",
//             option("--gfa").set(cfg.mode, output_type::gfa) % "produce graph in GFA1 format",
//             option("--spades").set(cfg.mode, output_type::spades) % "produce graph in SPAdes internal format",
//             option("--spades-gp").set(cfg.mode, output_type::spades_pack) % "produce graph pack in SPAdes internal format "
//                                                                        "(recommended if bulges are removed to improve further read mapping)")
  );

  auto result = parse(argc, argv, cli);
  if (!result) {
      std::cout << make_man_page(cli, argv[0]);
      exit(1);
  }
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

namespace debruijn {

namespace simplification {

class PositionStorage : public omnigraph::GraphActionHandler<Graph> {
    typedef Graph::EdgeId EdgeId;
    typedef Graph::VertexId VertexId;

    //FIXME support pos values other than -1
    std::multimap<EdgeId, int> poss_;

    void Insert(EdgeId e, int pos = -1) {
        poss_.insert(std::make_pair(e, pos));
    }

public:
    PositionStorage(const Graph &g) :
            omnigraph::GraphActionHandler<Graph>(g, "PositionStorage") {}

    void HandleDelete(EdgeId e) override {
        poss_.erase(e);
    }

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
        size_t cumm_len = 0;
        for (EdgeId e : old_edges) {
            for (int p : utils::get_all(poss_, e)) {
                Insert(new_edge, (p < 0) ? -1 : int(cumm_len + p));
            }
            cumm_len += g().length(e);
        }
    }

    void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) override {
        VERIFY_MSG(false, "No support");
    }

    void HandleSplit(EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) override {
        if (poss_.count(old_edge) == 0)
            return;
        if (old_edge == g().conjugate(old_edge)) {
            Insert(new_edge1);
            Insert(g().conjugate(new_edge1));
            Insert(new_edge2);
        } else {
            Insert(new_edge1);
            Insert(new_edge2);
        }
    }

    //FIXME support namer (in particular IdMapper)
    void Save(std::ostream &os) const {
        io::CanonicalEdgeHelper<Graph> canonical_helper(g());
        for (const auto &e_p : poss_) {
            //FIXME correct coord
            os << canonical_helper.EdgeOrientationString(e_p.first, "\t")
               << '\t' << e_p.second << '\n';
        }
    }

    void Load(std::istream &is, const io::EdgeLabelHelper<Graph> &label_helper) {
        std::string s;
        while (std::getline(is, s)) {
            std::istringstream ss(s);
            std::string label;
            ss >> label;
            EdgeId e = label_helper.edge(label);
            VERIFY_MSG(e != EdgeId(), "Couldn't find edge with int id " << std::stoi(label) << " in the graph");

            std::string orient;
            ss >> orient;
            VERIFY_MSG(orient == "+" || orient == "-", "Invalid orientation");

            //Currently ignored
            int pos;
            ss >> pos;

            e = (orient == "+") ? e : g().conjugate(e);
            Insert(e);
        }
    }

private:
    DECL_LOGGER("PositionStorage");
};

template<class Graph>
void FillFlankingFromAverage(const Graph &g, FlankingCoverage<Graph> &flanking_cov) {
    for (auto it = g.ConstEdgeBegin(); !it.IsEnd(); ++it) {
        EdgeId e = *it;
        auto raw_flank = g.coverage(e) * double(std::min(g.length(e), flanking_cov.averaging_range()));
        flanking_cov.SetRawCoverage(e, unsigned(math::round(raw_flank)));
    }
}

template<class Graph>
AlgoPtr<Graph> ConditionedTipClipperInstance(Graph &g,
                                             const config::debruijn_config::simplification::tip_clipper &tc_config,
                                             const SimplifInfoContainer &info,
                                             func::TypedPredicate<EdgeId> extra_condition,
                                             EdgeRemovalHandlerF<Graph> removal_handler = nullptr) {
    if (tc_config.condition.empty())
        return nullptr;

    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = func::And(parser(), extra_condition);
    auto algo = TipClipperInstance(g, condition, info, removal_handler);
    VERIFY_MSG(parser.requested_iterations() != 0, "To disable tip clipper pass empty string");
    if (parser.requested_iterations() == 1) {
        return algo;
    } else {
        return std::make_shared<LoopedAlgorithm<Graph>>(g, algo, 1, size_t(parser.requested_iterations()),
                /*force primary for all*/ true);
    }
}

static void Simplify(conj_graph_pack &gp,
                     const debruijn::simplification::SimplifInfoContainer &simplif_info,
                     bool rel_cov_proc_enabled) {

    const bool using_edge_qual = gp.edge_qual.IsAttached();

    const std::function<void(EdgeId)>& removal_handler = nullptr;

    typename ComponentRemover<Graph>::HandlerF set_removal_handler_f;
    if (removal_handler) {
        set_removal_handler_f = [=](const std::set<EdgeId> &edges) {
            std::for_each(edges.begin(), edges.end(), removal_handler);
        };
    }

    //Refill flanking coverage to get same behavior while working with gfa graphs
    FillFlankingFromAverage<Graph>(gp.g, gp.flanking_cov);
    VERIFY(gp.flanking_cov.IsAttached());

    INFO("Graph simplification started");
    size_t iteration = 0;
    auto message_callback = [&] () {
        INFO("PROCEDURE == Simplification cycle, iteration " << ++iteration);
    };

    config::debruijn_config::simplification::tip_clipper tc_config;
    tc_config.condition = "{ tc_lb 3.5, cb 1000000, rctc 2.0 }";

    omnigraph::CompositeAlgorithm<Graph> algo(gp.g, message_callback);

    func::TypedPredicate<EdgeId> extra_condition = func::AlwaysTrue<EdgeId>();
    if (gp.edge_qual.IsAttached()) {
        //FIXME do something to not remove undead ends.
        //FIXME Only in TipClipper or everywhere?
        extra_condition = [&] (EdgeId e) {
            return gp.edge_qual.IsZeroQuality(e);
        };
    }
    algo.AddAlgo(ConditionedTipClipperInstance(gp.g, tc_config, simplif_info,
                                               extra_condition, removal_handler),
                 "Tip clipper");

    //FIXME integrate condition
    algo.AddAlgo(BRInstance(gp.g, default_br_config(), simplif_info, removal_handler),
                        "Bulge remover");

    config::debruijn_config::simplification::erroneous_connections_remover ec_config;
    //ec_config.condition = "{ to_ec_lb 2, icb 2.5 }";
    ec_config.condition = "{ to_ec_lb 2, icb auto }";

    algo.AddAlgo(ECRemoverInstance(gp.g, ec_config, simplif_info, removal_handler),
                 "Low coverage edge remover with bounded length");

    AlgorithmRunningHelper<Graph>::IterativeThresholdsRun(algo,
                                                          /*cycle_iter_count*/3,
                                                          /*all_primary*/false);

    iteration = 0;

    auto low_cov_thr = std::max(2.0, simplif_info.detected_mean_coverage() / 100.);
    INFO("Unconditional coverage lower-bound set at " << low_cov_thr);
    //NB: we do not rescue the undeadends here
    algo.AddAlgo(std::make_shared<ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>>>
                        (gp.g,
                        CoverageUpperBound<Graph>(gp.g, low_cov_thr),
                        simplif_info.chunk_cnt(),
                        removal_handler,
                        /*canonical_only*/true,
                        CoverageComparator<Graph>(gp.g)),
                        "Removing all edges with coverage below " + std::to_string(low_cov_thr));

    algo.AddAlgo(RelativeCoverageComponentRemoverInstance<Graph>(gp.g, gp.flanking_cov,
            default_rcc_config(rel_cov_proc_enabled),
            simplif_info, set_removal_handler_f),
            "Removing subgraphs based on relative coverage");

    algo.AddAlgo(RelativelyLowCoverageDisconnectorInstance<Graph>(gp.g, gp.flanking_cov,
            default_red_config(rel_cov_proc_enabled),
            simplif_info, removal_handler),
            "Disconnecting relatively low covered edges");

    AlgorithmRunningHelper<Graph>::LoopedRun(algo, /*min it count*/1, /*max it count*/10);
    VERIFY(!using_edge_qual || gp.edge_qual.IsAttached());
    VERIFY(gp.flanking_cov.IsAttached());
}

}

}

//FIXME reduce code duplication with subgraph-extractor
static bool IsDeadEnd(const Graph &g, VertexId v) {
    return g.IncomingEdgeCount(v) * g.OutgoingEdgeCount(v) == 0;
}

//FIXME think about self-conjugate sources/sinks and other connections between RC subgrahps
//TODO improve, check consistency beteween total sink/source coverage estimates and think about weird cases
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

size_t DetermineSampleCnt(const std::string &profile_fn) {
    std::ifstream is(profile_fn);
    std::string line;
    std::getline(is, line);
    VERIFY_MSG(is, "I/O problem while reading " << profile_fn);
    std::istringstream ss(line);
    std::string token;
    size_t i = 0;
    while (ss >> token) {
        ++i;
    }
    VERIFY(i > 0);
    return i - 1;
}

//TODO set up reasonable flanking range
int main(int argc, char** argv) {
    utils::segfault_handler sh;
    gcfg cfg;

    process_cmdline(argc, argv, cfg);

    toolchain::create_console_logger();
    START_BANNER("SPAdes standalone graph simplifier");

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

        conj_graph_pack gp(k, tmpdir, 0);

        INFO("Loading de Bruijn graph from " << cfg.graph);
        omnigraph::GraphElementFinder<Graph> element_finder(gp.g);
        gp.kmer_mapper.Attach();

        io::EdgeLabelHelper<Graph> label_helper(element_finder,
                                                toolchain::LoadGraph(gp, cfg.graph));

        const Graph &g = gp.g;

        std::set<EdgeId> undeadends;
        if (!cfg.deadends_fn.empty()) {
            auto deadend_names = ReadDeadendNames(cfg.deadends_fn);
            for (auto it = g.ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
                EdgeId e = *it;
                if (IsDeadEnd(g, g.EdgeStart(e)) || IsDeadEnd(g, g.EdgeEnd(e))) {
                    if (!deadend_names.count(label_helper.label(e))) {
                        undeadends.insert(e);
                        undeadends.insert(g.conjugate(e));
                    }
                }
            }
            for (EdgeId e : undeadends) {
                VERIFY(undeadends.count(g.conjugate(e)));
                gp.edge_qual.AddQuality(e, 1000.);
            }
            gp.edge_qual.Attach();
        }

        using namespace debruijn_graph::coverage_profiles;
        typedef EdgeProfileStorage ProfileStorage;
        std::unique_ptr<ProfileStorage> profile_storage;

        if (!cfg.edge_profile_fn.empty()) {
            INFO("Loading edge profiles from " << cfg.edge_profile_fn);
            fs::CheckFileExistenceFATAL(cfg.edge_profile_fn);
            size_t sample_cnt = DetermineSampleCnt(cfg.edge_profile_fn);
            INFO("Sample count determined as " << sample_cnt);
            profile_storage = std::make_unique<ProfileStorage>(gp.g, sample_cnt);
            std::ifstream is(cfg.edge_profile_fn);
            profile_storage->Load(is, label_helper);
            INFO("Profiles loaded");
        }

        using debruijn::simplification::PositionStorage;
        std::unique_ptr<PositionStorage> stop_codons_storage;
        if (!cfg.stop_codons_fn.empty()) {
            INFO("Loading stop codon positions from " << cfg.stop_codons_fn);
            fs::CheckFileExistenceFATAL(cfg.stop_codons_fn);
            stop_codons_storage = std::make_unique<PositionStorage>(gp.g);
            std::ifstream is(cfg.stop_codons_fn);
            stop_codons_storage->Load(is, label_helper);
            INFO("Stop codon positions loaded");
        }

        debruijn::simplification::SimplifInfoContainer simplif_info;
        simplif_info.set_main_iteration(true);
        simplif_info.set_read_length(cfg.RL);
        simplif_info.set_chunk_cnt(200); //TODO
        double bin_cov = 0.;
        if (cfg.bin_cov_str.empty()) {
            INFO("Mean coverage option was not specified. Using default value of " << bin_cov);
        } else {
            INFO("Estimated mean coverage was specified as " << cfg.bin_cov_str)
            if (cfg.bin_cov_str == "auto") {
                INFO("Trying to determine from edge coverage");
                VERIFY_MSG(!cfg.deadends_fn.empty(),
                           "Deadends option (-d/--dead-ends) was not specified. "
                           "Can only determine coverage while working with subgraphs!");
                bin_cov = DetermineAvgCoverage(g, undeadends);
            } else {
                bin_cov = std::stod(cfg.bin_cov_str);
            }

        }
        simplif_info.set_detected_mean_coverage(bin_cov);
        auto ec_bound = std::max(2.5, bin_cov / 50.);
        simplif_info.set_detected_coverage_bound(ec_bound); //TODO
        INFO("Erroneous connection coverage bound set at " << ec_bound);

        debruijn::simplification::Simplify(gp, simplif_info, cfg.rel_cov_proc_enabled);

        INFO("Saving graph to " << cfg.outfile);

        fs::make_dirs(fs::parent_path(cfg.outfile));
        io::binary::BasePackIO<Graph>().Save(cfg.outfile, gp);

        if (cfg.save_gfa) {
            INFO("Saving gfa");
            std::ofstream os(cfg.outfile + ".gfa");
            gfa::GFAWriter writer(gp.g, os);
            writer.WriteSegmentsAndLinks();
        }

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
            VERIFY_MSG(gp.edge_qual.IsAttached(), "Edge quality got detached");
            std::ofstream deadends_os(cfg.outfile + ".deadends");
            for (auto it = g.ConstEdgeBegin(/*canonical_only*/true); !it.IsEnd(); ++it) {
                EdgeId e = *it;
                if (IsDeadEnd(g, g.EdgeStart(e)) || IsDeadEnd(g, g.EdgeEnd(e))) {
                    if (gp.edge_qual.IsPositiveQuality(e))
                        continue;
                    deadends_os << gp.g.int_id(e) << "\n";
                }
            }
        }

    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
