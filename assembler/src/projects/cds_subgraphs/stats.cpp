//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"
#include "toolchain/subgraph_utils.hpp"
#include "io/reads/file_reader.hpp"
#include "io/binary/graph_pack.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "utils/segfault_handler.hpp"

#include "version.hpp"

#include <clipp/clipp.h>
#include <string>
#include <set>
#include <numeric>
#include <sys/types.h>

using namespace debruijn_graph;

namespace debruijn_graph {

static void MappingPathStats(const Graph &g, const std::string &name, size_t query_kmer_len,
                             const omnigraph::MappingPath<EdgeId> &mapping_path,
                             io::EdgeNamingF<Graph> edge_naming_f) {
    VERIFY(mapping_path.size() > 0);
    INFO("Mapping path stats");
    std::stringstream ss;
    std::string delim = "";
    ss << "MappingPath ( ";
    for (size_t i = 0; i < mapping_path.size(); i++) {
        ss << delim << edge_naming_f(g, mapping_path[i].first) << ": " << mapping_path[i].second;
        delim = " ; ";
    }
    ss << " )";
    INFO(ss.str());

    size_t curr_pos = 0;
    for (size_t i = 0; i < mapping_path.size(); i++) {
        Range q_range = mapping_path[i].second.initial_range;
        //Range e_range = mapping_path[i].second.mapped_range;
        if (q_range.start_pos > curr_pos) {
            INFO("GAP while mapping " << name << " of length " << (q_range.start_pos - curr_pos) << " between " << curr_pos << " and " << q_range.start_pos);
        }
        curr_pos = q_range.end_pos;
    }
    if (query_kmer_len > curr_pos) {
        INFO("GAP while mapping " << name << " of length " << (query_kmer_len - curr_pos) << " between " << curr_pos << " and " << query_kmer_len);
    }
}

static bool UsedOnBothSides(const Graph &g, EdgeId e, const std::multimap<EdgeId, std::string> &edge_usage) {
    auto used_f = [&] (EdgeId a) { return edge_usage.count(a) > 0 || edge_usage.count(g.conjugate(a)) > 0; };
    return std::any_of(g.out_begin(g.EdgeEnd(e)), g.out_end(g.EdgeEnd(e)), used_f) &&
            std::any_of(g.in_begin(g.EdgeStart(e)), g.in_end(g.EdgeStart(e)), used_f);
}

static void Run(const GraphPack &gp, const std::string &contigs_file,
                const std::string &paths_fn, const std::string &edge_info_fn,
                const std::string &subgraph_prefix, const std::string &edge_color_fn,
                const io::EdgeLabelHelper<Graph> &label_helper) {

    const auto &graph = gp.get<Graph>();
    io::CanonicalEdgeHelper<Graph> canonical_helper(graph, label_helper.edge_naming_f());
    auto mapper = MapperInstance(gp);
    ReadPathFinder<Graph> path_finder(graph, /*skip_unfixed*/false);

    std::ofstream os;
    if (!paths_fn.empty())
        os.open(paths_fn);

    std::multimap<EdgeId, std::string> edge_usage;

    io::FileReadStream reader(contigs_file);
    io::SingleRead read;
    while (!reader.eof()) {
        reader >> read;
        if (!read.IsValid()) {
            WARN("Non-ACGT symbols in the " << read.name() << ". Not computing stats");
            continue;
        }

        INFO("Aligning " << read.name());
        os << read.name() << "\n";

        auto mapping_path = mapper->MapRead(read);
        if (mapping_path.empty()) {
            WARN("Failed to align " << read.name());
            continue;
        }

        MappingPathStats(graph, read.name(), read.sequence().size() - gp.k(), mapping_path, label_helper.edge_naming_f());
        std::string delimeter = "";
        auto fixed_path = path_finder.FindReadPath(mapping_path);
        if (!CheckContiguous(graph, fixed_path)) {
            INFO("Were not able to recover continuous path for " << read.name() << " even after fixing attempt");
        }
        for (const auto& e : fixed_path) {
            edge_usage.insert(make_pair(canonical_helper.Canonical(e), read.name()));
            os << delimeter << canonical_helper.EdgeOrientationString(e);
            delimeter = ",";
        }
        os << "\n";
    }

    std::stringstream suspicious_edges;
    std::ofstream edge_info_os;
    if (!edge_info_fn.empty())
        edge_info_os.open(edge_info_fn);

    std::ofstream edge_color_os;
    if (!edge_color_fn.empty()) {
        edge_color_os.open(edge_color_fn);
        edge_color_os << "Name,color\n";
    }

    //TODO optimize if we do not need to output unused edges
    if (!edge_info_fn.empty() || !edge_color_fn.empty()) {
        size_t suspicious_edge_cnt = 0;
        size_t total_suspicious_kmer_len = 0;

        for (EdgeId e : graph.canonical_edges()) {
            edge_info_os << label_helper.label(e);
            std::string delimeter = "\t";
            if (edge_usage.count(e) == 0 && UsedOnBothSides(graph, e, edge_usage)) {
                suspicious_edge_cnt++;
                suspicious_edges << label_helper.label(e) << ",";
                total_suspicious_kmer_len += graph.length(e);
                edge_color_os << label_helper.label(e) << ",red\n";
            }
            if (edge_usage.count(e) > 0) {
                edge_color_os << label_helper.label(e) << ",green\n";
            }

            for (const auto &usage : utils::get_all(edge_usage, e)) {
                edge_info_os << delimeter << usage;
            }
            edge_info_os << "\n";
        }

        INFO("Suspicious edge cnt: " << suspicious_edge_cnt << "; " <<
                "total (" << (gp.k() + 1) << "-mer) length: " << total_suspicious_kmer_len);
        INFO("Suspicious edges: " << suspicious_edges.str());
    }

    if (!subgraph_prefix.empty()) {
        auto edges = utils::key_set(edge_usage);
        toolchain::ComponentExpander expander(graph);
        auto component = expander.Expand(omnigraph::GraphComponent<Graph>::FromEdges(graph, edges.begin(),
                                                                          edges.end(), /*add conjugate*/true));
        toolchain::WriteComponentWithDeadends(component, subgraph_prefix, label_helper.edge_naming_f());
    }
}

} // namespace debruijn_graph

struct gcfg {
    gcfg() : k(0),
             workdir("./tmp"),
             nthreads(omp_get_max_threads() / 2 + 1)
    {}

    unsigned k;
    std::string graph_path;
    std::string sequences_fn;
    std::string paths_fn;
    std::string edge_info_fn;
    std::string edge_color_fn;
    std::string subgraph_prefix;
    std::string workdir;
    unsigned nthreads;
};

static void process_cmdline(int argc, char **argv, gcfg &cfg) {
    using namespace clipp;

    auto cli = (
               (required("-k") & integer("value", cfg.k)) % "k-mer length to use",
               (required("-g", "--graph") & value("graph", cfg.graph_path)) % "In GFA (ending with .gfa) or prefix to SPAdes graph pack",
               (required("-q", "--cds-queries") & value("file", cfg.sequences_fn)) % "Path to FASTA file with ground truth CDS sequences",
               (option("-p", "--paths") & value("file", cfg.paths_fn)) % "Destination for outputting paths corresponding to CDS sequences",
               (option("-e", "--edge-info") & value("file", cfg.edge_info_fn)) % "Destination for outputting edge usage information",
               (option("-t", "--threads") & integer("value", cfg.nthreads)) % "# of threads to use (default: max_threads / 2)",
               (option("-c", "--colors") & value("file", cfg.edge_color_fn)) % "Destination for outputting edge coloring to be displayed in Bandage",
               (option("-s", "--subgraph") & value("file", cfg.subgraph_prefix)) % "Destination for outputting locality of covered edges in GFA",
               (option("--workdir") & value("dir", cfg.workdir)) % "scratch directory to use (default: ./tmp)"
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
    START_BANNER("Analyzing alignment of known CDS sequences to assembly graph or its subgraph");

    try {
        unsigned k = cfg.k;
        INFO("K-mer length set to " << k);

        unsigned nthreads = cfg.nthreads;
        nthreads = std::min(nthreads, (unsigned) omp_get_max_threads());
        // Inform OpenMP runtime about this :)
        omp_set_num_threads((int) nthreads);
        INFO("# of threads to use: " << nthreads);

        fs::make_dirs(cfg.workdir);
        debruijn_graph::GraphPack gp(k, cfg.workdir, 0);

        omnigraph::GraphElementFinder<Graph> element_finder(gp.get<Graph>());
        INFO("Loading de Bruijn graph from " << cfg.graph_path);
        gp.get_mutable<debruijn_graph::KmerMapper<debruijn_graph::Graph>>().Attach(); // TODO unnecessary
        io::EdgeLabelHelper<Graph> label_helper(element_finder,
                toolchain::LoadGraph(gp, cfg.graph_path));

        gp.EnsureBasicMapping();

        debruijn_graph::Run(gp, cfg.sequences_fn, cfg.paths_fn, cfg.edge_info_fn,
                            cfg.subgraph_prefix, cfg.edge_color_fn, label_helper);
        INFO("Done");
    } catch (const std::string &s) {
        std::cerr << s << std::endl;
        return EINTR;
    } catch (const std::exception &e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return EINTR;
    }
}
