#include "stats.hpp"

#include "common/io/graph/gfa_reader.hpp"
#include "common/utils/logger/log_writers.hpp"

#include <filesystem>
#include <string>
#include <unordered_map>

#include <clipp/clipp.h>


void ReadGraph(debruijn_graph::Graph& g, const std::filesystem::path& read_from) {
    gfa::GFAReader reader(read_from);
    reader.to_graph(g);
}

void CreateConsoleLogger(const std::filesystem::path& log_fn="") {
    using namespace logging;
    logger *lg = create_logger(exists(log_fn) ? log_fn : "");
    lg->add_writer(std::make_shared<console_writer>());
    attach_logger(lg);
}

int main(int argc, char *argv[]) {
    CreateConsoleLogger();
    std::string command;
    std::string graph_path_;
    std::vector<double> n_percentiles_;
    std::vector<double> median_length_percentiles_;
    std::vector<double> sequences_coverage_percentiles_;

    using namespace clipp;

    auto cli = (
            value("command name", command),
            value("path to the graph", graph_path_),
            option("-n", "--n_percentiles") & numbers("Percentiles for Nx stat", n_percentiles_),
            option("-l", "--median_length_percentiles") & numbers("Median length for a given percentile of sequences", median_length_percentiles_),
            option("-v", "--sequences_coverage_percentiles") & numbers("The percentiles for depth of the graph, by base", sequences_coverage_percentiles_)
    );
    parse(argc, argv, cli);
/*
    std::vector<double> n_percentiles;
    for (auto perc_ : n_percentiles_) {
        n_percentiles.push_back(std::stod(perc_));
        std::cout << perc_ << "\n";
    }
    std::vector<double> median_length_percentiles;
    for (auto perc_ : median_length_percentiles_) {
        median_length_percentiles.push_back(std::stod(perc_));
        std::cout << perc_ << "\n";
    }
    std::vector<double> sequences_coverage_percentiles;
    for (auto perc_ : sequences_coverage_percentiles_) {
        sequences_coverage_percentiles.push_back(std::stod(perc_));
        std::cout << perc_ << "\n";
    }
*/
    std::filesystem::path graph_path = graph_path_;
    debruijn_graph::Graph g(55); //это ок, что я не знаю k, но создать граф без него не могу?
    ReadGraph(g, graph_path);

    if (command == "stat") {
        PrintStat(g, n_percentiles_, median_length_percentiles_, sequences_coverage_percentiles_);
    }

    return 0;
}