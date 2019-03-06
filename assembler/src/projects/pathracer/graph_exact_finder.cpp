//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "assembly_graph/core/graph.hpp"
#include "fasta_reader.hpp"
#include "io/graph/gfa_reader.hpp"
#include "naive_edge_index.hpp"

#include <fstream>
#include "common/utils/verify.hpp"
#include "io/reads/osequencestream.hpp"

#include <clipp/clipp.h>

#include "utils/logger/log_writers.hpp"
#include "utils/logger/log_writers_thread.hpp"
#include "utils.hpp"

void create_console_logger(const std::string &filename = "") {
    using namespace logging;

    logger *lg = create_logger("");
    lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<console_writer_thread>()));
    if (filename != "") {
        lg->add_writer(std::make_shared<mutex_writer>(std::make_shared<file_writer_thread>(filename)));
    }
    attach_logger(lg);
}

enum class Mode {
    none,
    nt,
    aa
};

int main(int argc, char* argv[]) {
    std::string graph_file;
    std::string sequence_file;
    std::string output_file;
    size_t k;
    Mode mode;
    using namespace clipp;
    auto cli =
        (sequence_file << value("input sequence file"),
         graph_file << value("graph file in GFA"),
         k << integer("k-mer size"),
         required("--output", "-o") & value("output file", output_file) % "output file",
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
    INFO("Edge index construction");
    NaiveEdgeIndex eindex;
    if (mode == Mode::aa) {
        eindex.build_aa(graph, (k - 2) / 3);
    } else {
        eindex.build_nt(graph, k);
    }
    INFO("Index built");

    std::ofstream of(output_file);
    for (auto &record : read_fasta(sequence_file)) {
        auto &seq = record.second;
        if (seq.find('X') != std::string::npos) {
            INFO("Stop-codon found, skipping the sequence");
            continue;
        }
        // remove -
        seq.erase(std::remove(seq.begin(), seq.end(), '-'), seq.end());
        to_upper_case(seq);
        INFO("Checking >" << record.first);
        INFO(seq);
        size_t count = mode == Mode::aa ? eindex.has_aa(seq, &graph) : eindex.has_nt(seq, &graph);
        INFO("#" << eindex.has_aa(seq, &graph));
        if (count) {
            of << ">" << record.first << "|COUNT=" << count <<  "\n";
            io::WriteWrapped(seq, of);
        }
    }

    return 0;
}

