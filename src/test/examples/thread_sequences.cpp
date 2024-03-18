//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alignment/sequence_mapper.hpp"
#include "io/reads/file_reader.hpp"
#include "io/binary/graph_pack.hpp"
#include "pipeline/graph_pack_helpers.h"
#include "pipeline/sequence_mapper_gp_api.hpp"
#include "toolchain/edge_label_helper.hpp"
#include "toolchain/utils.hpp"
#include "utils/segfault_handler.hpp"

#include <clipp/clipp.h>

using namespace std;

namespace debruijn_graph {

static void Run(size_t K, const filesystem::path &graph_path, const filesystem::path &contigs_file,
         const filesystem::path& out_paths_fn, const filesystem::path& out_edge_info_fn,
         const filesystem::path &tmpdir) {
    create_directory(tmpdir);

    graph_pack::GraphPack gp(K, tmpdir, 0);
    const auto& graph = gp.get<Graph>();
    INFO("Loading de Bruijn graph from " << graph_path);
    omnigraph::GraphElementFinder<Graph> element_finder(graph);
    gp.get_mutable<KmerMapper<Graph>>().Attach();
    io::EdgeLabelHelper<Graph> label_helper(element_finder,
                                            toolchain::LoadBaseGraph(gp, graph_path));

    EnsureBasicMapping(gp);

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

    using namespace clipp;

    try {
        unsigned k = 55;
        std::filesystem::path graph_path, sequences_fn, path_out_fn, edge_out_fn;
        std::string graph, sequences, path_out, edge_out, workdir = "tmp";
        
        auto cli = (
            required("-k", "--kmer") & value("k", k) % "K-mer length",
            required("-g", "--graph") & value("graph", graph) % "GFA file or folder with SPAdes saves",
            required("-q" "--queries") & value("file", sequences) % "Query sequences",
            required("-p", "--path_out_file") & value("paths", path_out) % "File to store path output",
            required("-e", "--edge_out_file") & value("edges", edge_out) % "File to store edge output",
            option("-w", "--workdir") & value("workdir", workdir) % "Wording directory (default: ./tmp)"
        );

        if (!parse(argc, argv, cli)) {
            std::cout << make_man_page(cli, argv[0]);
            exit(1);
        }
        
        graph_path = graph; path_out_fn = path_out; edge_out_fn = edge_out;

        toolchain::create_console_logger();

        START_BANNER("Threading sequences through assembly graph");

        INFO("K-mer length set to " << k);

        debruijn_graph::Run(k, graph_path, sequences_fn, path_out_fn, edge_out_fn, workdir);
    } catch (const std::string &s) {
        std::cerr << s;
        return EINTR;
    }

    return 0;
}
