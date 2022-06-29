#pragma once

#include <filesystem>
#include <optional>
#include <string>
#include <vector>


class CommandLineArguments {
 public:
    CommandLineArguments(int argc, char *argv[]);
    ~CommandLineArguments() = default;

    const std::string GetCommand() const { return command;}
    const std::filesystem::path GetGraphPath() const {return graph_path;}
    const std::vector<double> GetNPercentiles() const {return n_percentiles;}
    const std::vector<double> GetLengthPercentile() const {return median_length_percentiles; }
    const std::vector<double> GetSequencesCoveragePercentiles() const {return sequences_coverage_percentiles;}
    const std::filesystem::path GetYamlOutputPath() const {return yaml_output_path; }
    const std::optional<double> GetMinimalCoverageThreshold() const {return minimal_coverage_threshold;}
    const std::optional<double> GetMaximalCoverageThreshold() const {return maximal_coverage_threshold;}
    const std::optional<size_t> GetMinimalLengthThreshold() const {return minimal_length_threshold;}
    const std::optional<size_t> GetMaximalLengthThreshold() const {return maximal_length_threshold;}
    const bool GetDeleteIsolatedEdges() const {return delete_isolated_edges;}
    const std::filesystem::path GetGraphOutputPath() const {return graph_output_path;}

 private:
    std::string command;
    std::filesystem::path graph_path;

    std::vector<double> n_percentiles;
    std::vector<double> median_length_percentiles;
    std::vector<double> sequences_coverage_percentiles;
    std::filesystem::path yaml_output_path;

    std::optional<double> maximal_coverage_threshold;
    std::optional<double> minimal_coverage_threshold;
    std::optional<size_t> maximal_length_threshold;
    std::optional<size_t> minimal_length_threshold;
    bool delete_isolated_edges = false;

    std::filesystem::path graph_output_path;
};
