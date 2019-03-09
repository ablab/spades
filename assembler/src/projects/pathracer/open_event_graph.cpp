//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "io/graph/gfa_reader.hpp"

#include <fstream>
#include "common/utils/verify.hpp"

#include <clipp/clipp.h>

#include "utils.hpp"
#include "utils/logger/log_writers.hpp"
#include "utils/logger/log_writers_thread.hpp"

#include <cereal/archives/binary.hpp>
#include <fstream>

#include "pathtree.hpp"
#include "aa_cursor.hpp"
#include "cached_cursor.hpp"
#include "omnigraph_wrapper.hpp"

void create_console_logger(const std::string &filename = "") {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer_thread>()));
    if (filename != "") {
        lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<file_writer_thread>(filename)));
    }
    attach_logger(lg);
}

enum class Mode { none, nt, aa };

int main(int argc, char *argv[]) {
    std::string graph_file;
    std::string cereal_file;
    // std::string sequence_file;
    // std::string output_file;
    size_t k;
    Mode mode;
    using namespace clipp;
    auto cli =
        (
         cereal_file << value("cerealized event graph"),
         graph_file << value("graph file in GFA"),
         k << integer("k-mer size"),
         // required("--output", "-o") & value("output file", output_file) % "output file",
         (option("--nt").set(mode, Mode::nt) % "match against nucleotide string(s)" |
          option("--aa").set(mode, Mode::aa) % "match agains amino acid string(s)"));

    if (!parse(argc, argv, cli)) {
        std::cout << make_man_page(cli, argv[0]);
        exit(1);
    }

    create_console_logger();
    VERIFY(mode != Mode::none);

    INFO("Reading GFA graph");
    debruijn_graph::ConjugateDeBruijnGraph graph(k);
    gfa::GFAReader gfa(graph_file);
    INFO("GFA segments: " << gfa.num_edges() << ", links: " << gfa.num_links());
    gfa.to_graph(graph, nullptr);
    INFO("Graph loaded");

    size_t top = 100;

    if (mode == Mode::aa) {
        std::vector<AAGraphCursor<DebruijnGraphCursor>> cursors;
        CachedCursorContext ccc(cursors, &graph);
        PathSet<CachedCursor> result(nullptr);
        std::ifstream ifs(cereal_file);
        cereal::BinaryInputArchive iarchive(ifs);
        iarchive(cursors, ccc, result);

        INFO("Extracting top paths");
        auto top_paths = result.top_k(top);
        bool x_as_m_in_alignment = mode == Mode::aa;
        if (!top_paths.empty()) {
            INFO("Best score in the current component: " << result.best_score());
            INFO("Best sequence in the current component");
            INFO(top_paths.str(0, &ccc));
            INFO(top_paths.str(1, &ccc));
            // INFO("Alignment: " << compress_alignment(top_paths.alignment(0, fees, &ccc), x_as_m_in_alignment));
        }
    }

    return 0;
}

