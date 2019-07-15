#include "transition_subgraph_validation.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

boost::optional<vector<SimpleTransitionGraphValidator::ScaffoldVertex>> SimpleTransitionGraphValidator::GetCorrectPath(
    const SimpleTransitionGraphValidator::SimpleTransitionGraph &graph,
    scaffold_graph::ScaffoldVertex source,
    scaffold_graph::ScaffoldVertex sink) const {
    DEBUG("Getting correct path");
    boost::optional<vector<ScaffoldVertex>> result;
    vector<ScaffoldVertex> intermediate_result;
    ScaffoldVertex current = source;
    bool got_next = false;
    intermediate_result.push_back(current);
    set<ScaffoldVertex> visited;
    while (current != sink and visited.find(current) == visited.end()) {
        visited.insert(current);
        for (auto it = graph.outcoming_begin(current); it != graph.outcoming_end(current); ++it) {
            auto next = *it;
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
    GeneralTransitionStorageBuilder reference_storage_builder(gp_.g, 1, false, false);
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold_);
    auto reference_transition_storage = reference_storage_builder.GetTransitionStorage(reference_paths);
    SimpleTransitionGraphValidator transition_graph_validator(reference_transition_storage, length_threshold_);
    return transition_graph_validator;
}
SimpleTransitionGraphValidatorConstructor::SimpleTransitionGraphValidatorConstructor(const conj_graph_pack &gp,
                                                                                     size_t length_threshold) :
    gp_(gp), length_threshold_(length_threshold) {}
}
}
}