//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "transition_subgraph_validation.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

boost::optional<std::vector<scaffold_graph::ScaffoldVertex>> SimpleTransitionGraphValidator::GetCorrectPath(
    const SimpleTransitionGraphValidator::SimpleTransitionGraph &graph,
    scaffold_graph::ScaffoldVertex source,
    scaffold_graph::ScaffoldVertex sink) const {
    DEBUG("Getting correct path");
    boost::optional<std::vector<ScaffoldVertex>> result;
    std::vector<ScaffoldVertex> intermediate_result;
    ScaffoldVertex current = source;
    bool got_next = false;
    intermediate_result.push_back(current);
    std::set<ScaffoldVertex> visited;
    while (current != sink and visited.find(current) == visited.end()) {
        visited.insert(current);
        for (const auto &next: graph.OutNeighbours(current)) {
            VERIFY(graph.ContainsVertex(next));
            TRACE("Checking transition");
            bool is_correct = reference_transition_storage_.CheckTransition(current, next);
            TRACE("Checked transition");
            if (is_correct) {
                intermediate_result.push_back(next);
                TRACE("Current vertex: " << current.int_id());
                TRACE("Next vertex: " << next.int_id())
                current = next;
                got_next = true;
                break;
            }
        }
        if (not got_next) {
            return result;
        }
    }
    if (not(current == sink)) {
        return result;
    }
    result = intermediate_result;
    return result;
}
SimpleTransitionGraphValidator::SimpleTransitionGraphValidator(const ContigTransitionStorage &reference_transition_storage,
                                                               size_t length_threshold)
    : reference_transition_storage_(reference_transition_storage), length_threshold_(length_threshold) {}
SimpleTransitionGraphValidator SimpleTransitionGraphValidatorConstructor::GetValidator(
    const string &path_to_reference) const {
    GeneralTransitionStorageBuilder reference_storage_builder(g_, 1, false, false);
    validation::FilteredReferencePathHelper path_helper(g_, index_, kmer_mapper_);
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, length_threshold_, 5000.0);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromUnique(path_to_reference, unique_storage);
    auto reference_transition_storage = reference_storage_builder.GetTransitionStorage(reference_paths);
    SimpleTransitionGraphValidator transition_graph_validator(reference_transition_storage, length_threshold_);
    return transition_graph_validator;
}
SimpleTransitionGraphValidatorConstructor::SimpleTransitionGraphValidatorConstructor(
        const GraphPack &gp,
        size_t length_threshold) :
    gp_(gp),
    g_(gp_.get<Graph>()),
    index_(gp_.get<EdgeIndex<Graph>>()),
    kmer_mapper_(gp_.get<KmerMapper<Graph>>()),
    length_threshold_(length_threshold) {}
}
}
}