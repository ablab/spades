//***************************************************************************
//* Copyright (c) 2020 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "pipeline_validation.hpp"

#include "scaffold_graph_validation.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

void ScaffoldGraphPipelineValidator::ValidateStagesResults(const Pipeline &pipeline,
                                                           const std::string &output_path) const {
    const auto intermediate_results = pipeline.GetIntermediateResults();
//    const string path_to_reference = configs_.statistics.genome_path;
    DEBUG("Path to reference: " << path_to_reference_);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference_));
    const auto &graph = gp_.get<Graph>();
    const auto &index = gp_.get<EdgeIndex<Graph>>();
    const auto &kmer_mapper = gp_.get<KmerMapper<Graph>>();
    validation::ScaffoldGraphValidator scaffold_graph_validator(graph);
    validation::FilteredReferencePathHelper path_helper(graph, index, kmer_mapper);
    VERIFY_DEV(not intermediate_results.empty());
    auto reference_paths = path_helper.GetFilteredReferencePathsFromUnique(path_to_reference_,
                                                                           unique_storage_);

    for (const auto &result: intermediate_results) {
        auto scaffold_graph = *(result.first);
        string name = result.second;
        auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
        INFO("Stats for " << name);
        stats.Serialize(std::cout, false);
        std::string stage_path = fs::append_path(output_path, name);
        std::ofstream fout(stage_path);
        stats.Serialize(fout, true);
    }
}
ScaffoldGraphPipelineValidator::ScaffoldGraphPipelineValidator(const string &path_to_reference,
                                                               const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                               const GraphPack &gp) :
    path_to_reference_(path_to_reference),
    unique_storage_(unique_storage),
    gp_(gp) {}

}
}
}