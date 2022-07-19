#include "parser.hpp"

#include <clipp/clipp.h>

#include <cassert>

#include <iostream>

CommandLineArguments::CommandLineArguments(int argc, char *argv[]) {
    std::string graph_path_;
    std::string yaml_output_path_;
    std::string graph_output_path_;

    bool is_minimal_coverage_threshold = false;
    double minimal_coverage_threshold_ = 10;

    bool is_maximal_coverage_threshold;
    double maximal_coverage_threshold_ = 1000;

    bool is_maximal_length_threshold = false;
    size_t maximal_length_threshold_ = 10000;

    bool is_minimal_length_threshold = false;
    size_t minimal_length_threshold_ = 50;


    using namespace clipp;
    auto cli = (
            value("command name", command),
            value("path to the graph", graph_path_),
            option("-n", "--n_percentiles") & numbers("Percentiles for Nx stat", n_percentiles),
            option("-l", "--median_length_percentiles") & numbers("Median length for a given percentile of sequences", median_length_percentiles),
            option("-v", "--sequences_coverage_percentiles") & numbers("The percentiles for depth of the graph, by base", sequences_coverage_percentiles),
            option("-y", "--yaml_output") & value("YAML output file", yaml_output_path_),
            option("--mincoverage").set(is_minimal_coverage_threshold) & opt_value("minimal coverage", minimal_coverage_threshold_),
            option("--maxcoverage").set(is_maximal_coverage_threshold) & opt_value("maximal coverage", maximal_coverage_threshold_),
            option("--maxlen").set(is_maximal_length_threshold) & value("maximal sequence length", maximal_length_threshold_),
            option("--minlen").set(is_minimal_length_threshold) & value("minimal sequence length", minimal_length_threshold_),
            option("--del_isolated").set(delete_isolated_edges),
            option("--output_graph") & value("path to output graph", graph_output_path_),
            option("--sequence_scope") & value("name of a sequence to find a scope around", sequence_name_to_find_scope),
            option("--scope_depth") & value("depth of a scope", scope_depth)
    );
    parse(argc, argv, cli);
    graph_path = graph_path_;
    assert(exists(graph_path) && "File doesn't exist");
    yaml_output_path = yaml_output_path_;

    if (is_minimal_coverage_threshold) {
        minimal_coverage_threshold = minimal_coverage_threshold_;
    }
    if (is_maximal_coverage_threshold) {
        maximal_coverage_threshold = maximal_coverage_threshold_;
    }
    if (is_maximal_length_threshold) {
        maximal_length_threshold = maximal_length_threshold_;
    }
    if (is_minimal_length_threshold) {
        minimal_length_threshold = minimal_length_threshold_;
    }

    std::sort(n_percentiles.begin(), n_percentiles.end());
    std::sort(median_length_percentiles.begin(), median_length_percentiles.end());
    std::sort(sequences_coverage_percentiles.begin(), sequences_coverage_percentiles.end());

    graph_output_path = graph_output_path_;
}