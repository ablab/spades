//***************************************************************************
//* Copyright (c) 2018 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"
#include "io/reads/file_reader.hpp"
#include "io/binary/graph_pack.hpp"
#include "modules/alignment/sequence_mapper.hpp"
#include "utils/segfault_handler.hpp"

#include <cxxopts/cxxopts.hpp>

using namespace std;

namespace debruijn_graph {

static void Run(size_t K, const string &graph_path, const string &contigs_file,
         const string& out_paths_fn, const string& out_edge_info_fn,
         const string &tmpdir) {
    fs::make_dir(tmpdir);

    GraphPack gp(K, tmpdir, 0);
    const auto& graph = gp.get<Graph>();
    INFO("Loading de Bruijn graph from " << graph_path);
    omnigraph::GraphElementFinder<Graph> element_finder(graph);
    gp.get_mutable<KmerMapper<Graph>>().Attach();
    io::EdgeLabelHelper<Graph> label_helper(element_finder,
                                            toolchain::LoadGraph(gp, graph_path));

    gp.EnsureBasicMapping();

    ReadPathFinder<Graph> path_finder(graph, /*skip_unfixed*/false);
    auto mapper = MapperInstance(gp);

    io::FileReadStream reader(contigs_file);
    io::CanonicalEdgeHelper<Graph> canonical_helper(graph, label_helper.edge_naming_f());

    std::ofstream os(out_paths_fn);
    std::multimap<EdgeId, std::string> edge_usage;
    io::SingleRead read;
    while (!reader.eof()) {
        reader >> read;
        INFO("Aligning " << read.name());
        os << read.name() << "\n";

        std::string delimeter = "";
        for (const auto& e : path_finder.FindReadPath(mapper->MapRead(read))) {
            edge_usage.insert(make_pair(canonical_helper.Canonical(e), read.name()));
            os << delimeter << canonical_helper.EdgeOrientationString(e);
            delimeter = ",";
        }
        os << "\n";
    }

    std::ofstream edge_info_os(out_edge_info_fn);
    for (EdgeId e : graph.canonical_edges()) {
        edge_info_os << label_helper.label(e);
        std::string delimeter = "\t";
        for (const auto &usage : utils::get_all(edge_usage, e)) {
            edge_info_os << delimeter << usage;
        }
        edge_info_os << "\n";
    }
}

}

int main(int argc, char** argv) {
    srand(42);
    srandom(42);

    try {
        unsigned k;
        std::string workdir, graph_path, sequences_fn, path_out_fn, edge_out_fn;

        cxxopts::Options options(argv[0], " thread sequences through the graph");
        options.add_options()
                ("k,kmer", "K-mer length", cxxopts::value<unsigned>(k)->default_value("55"), "K")
                ("g,graph", "GFA file or folder with SPAdes saves", cxxopts::value<std::string>(graph_path))
                ("q,queries", "Query sequences", cxxopts::value<std::string>(sequences_fn), "file")
                ("p,path_out_file", "File to store path output", cxxopts::value<std::string>(path_out_fn))
                ("e,edge_out_file", "File to store edge output", cxxopts::value<std::string>(edge_out_fn))
                ("w,workdir", "Working directory (default: ./tmp)", cxxopts::value<std::string>(workdir)->default_value("./tmp"), "dir")
                ("h,help", "Print help");

        options.parse(argc, argv);
        if (options.count("help")) {
            std::cout << options.help() << std::endl;
            exit(0);
        }

        toolchain::create_console_logger();

        START_BANNER("Threading sequences through assembly graph");

        INFO("K-mer length set to " << k);

        debruijn_graph::Run(k, graph_path, sequences_fn, path_out_fn, edge_out_fn, workdir);
    } catch (const std::string &s) {
        std::cerr << s;
        return EINTR;
    } catch (const cxxopts::OptionException &e) {
        std::cerr << "error parsing options: " << e.what() << std::endl;
        exit(1);
    }
}
