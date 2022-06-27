//***************************************************************************
//* Copyright (c) 2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "stats.hpp"

#include "common/io/graph/gfa_reader.hpp"
#include "common/toolchain/utils.hpp"

#include <clipp/clipp.h>

int main(int argc, char *argv[]) {
    toolchain::create_console_logger();

    std::string command;
    std::string graph_path_;
    std::string yaml_output_path_;
    std::vector<double> n_percentiles;
    std::vector<double> median_length_percentiles;
    std::vector<double> sequences_coverage_percentiles;

    using namespace clipp;

    auto cli = (
            value("command name", command),
            value("path to the graph", graph_path_),
            option("-n", "--n_percentiles") & numbers("Percentiles for Nx stat", n_percentiles),
            option("-l", "--median_length_percentiles") & numbers("Median length for a given percentile of sequences", median_length_percentiles),
            option("-v", "--sequences_coverage_percentiles") & numbers("The percentiles for depth of the graph, by base", sequences_coverage_percentiles),
            option("-y", "--yaml_output") & value("YAML output file", yaml_output_path_)
    );
    parse(argc, argv, cli);

    std::filesystem::path graph_path = graph_path_;
    std::filesystem::path yaml_output_path = yaml_output_path_;

    gfa::GFAReader reader(graph_path);
    debruijn_graph::Graph g(reader.k());
    reader.to_graph(g);

    if (command == "stat") {
        std::sort(n_percentiles.begin(), n_percentiles.end());
        std::sort(median_length_percentiles.begin(), median_length_percentiles.end());
        std::sort(sequences_coverage_percentiles.begin(), sequences_coverage_percentiles.end());
        gfa_tools::CalculateStat(g, n_percentiles, median_length_percentiles, sequences_coverage_percentiles, yaml_output_path);
    }

    return 0;
}