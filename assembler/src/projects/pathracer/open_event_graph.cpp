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
#include "io/reads/osequencestream.hpp"

#include "pathtree.hpp"
#include "aa_cursor.hpp"
#include "cached_cursor.hpp"
#include "omnigraph_wrapper.hpp"

#include "cursor_conn_comps.hpp"  // For to_nucl_path --- FIXME factor out this fnc

void create_console_logger(const std::string &filename = "") {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer_thread>()));
    if (filename != "") {
        lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<file_writer_thread>(filename)));
    }
    attach_logger(lg);
}

// FIXME factor them out

template <typename... Ts> using void_t = void;

template <typename T, typename = void>
struct has_edge_method : std::false_type {};

template <typename T>
struct has_edge_method<T, void_t<decltype(T{}.edge())>> : std::true_type {};
template<class GraphCursor>
std::enable_if_t<has_edge_method<GraphCursor>::value, std::vector<typename GraphCursor::EdgeId>> to_path(const std::vector<GraphCursor> &cpath) {
    std::vector<typename GraphCursor::EdgeId> path;

    size_t prev_position = 0;
    for (auto cursor : cpath) {
        const auto e = cursor.edge();
        size_t position = cursor.position();
        if (path.empty() || e != path.back() || prev_position >= position) {
            path.push_back(e);
        }
        prev_position = position;
    }

    return path;
}

template<class GraphCursor>
std::vector<typename GraphCursor::EdgeId> to_path(const std::vector<AAGraphCursor<GraphCursor>> &cpath) {
    return to_path(to_nucl_path(cpath));
}

using debruijn_graph::EdgeId;
using debruijn_graph::VertexId;
using debruijn_graph::ConjugateDeBruijnGraph;
enum class Mode { none, nt, aa };

int main(int argc, char *argv[]) {
    std::string graph_file;
    std::string cereal_file;
    std::string output_file;
    size_t top = 100;
    // std::string hmm_file;
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
         required("--output", "-o") & value("output file", output_file)    % "output file",
         (option("--top") & integer("N", top)) % "extract top N paths",
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

    std::ofstream of(output_file);

    auto process = [&](const auto &cursor) {
        using Cursor = std::decay_t<decltype(cursor)>;
        std::vector<Cursor> cursors;
        CachedCursorContext ccc(cursors, &graph);
        PathSet<CachedCursor> result(nullptr);
        std::ifstream ifs(cereal_file);
        cereal::BinaryInputArchive iarchive(ifs);
        iarchive(cursors, ccc, result);

        // INFO("Collapsing");
        // const_cast<pathtree::PathLink<CachedCursor>*>(result.pathlink())->apply([](auto &pl) { pl.collapse_and_trim(); });

        INFO("Extracting top paths");
        auto top_paths = result.top_k(top);
        // bool x_as_m_in_alignment = mode == Mode::aa;
        VERIFY(!top_paths.empty());

        size_t count = 0;
        for (const auto& annotated_path : top_paths) {
            VERIFY(annotated_path.path.size());
            std::string seq = annotated_path.str(&ccc);
            auto unpacked_path = ccc.UnpackPath(annotated_path.path, cursors);
            auto nucl_path = to_nucl_path(unpacked_path);
            auto edge_path = to_path(nucl_path);
            if (count % 10000 == 0) {
                size_t pos = nucl_path[0].position();
                INFO(annotated_path.score);
                INFO(seq);
                INFO(edge_path << " position = " << pos);
            }
            ++count;
            of << ">Score=" << annotated_path.score << "|edges=" << edge_path << "\n";
            io::WriteWrapped(seq, of);
        }
    };
    if (mode == Mode::aa) {
        process(AAGraphCursor<DebruijnGraphCursor>());
    } else {
        process(DebruijnGraphCursor());
    }

    return 0;
}

