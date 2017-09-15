//
// Created by itolstoganov on 06.10.17.
//

#include "pe_extraction.hpp"

namespace path_extend {
bool SimplePEPredicate::Check(const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const {
    return transition_storage_.CheckTransition(scaffold_edge.getStart(), scaffold_edge.getEnd());
}
SimplePEPredicate SimplePEPredicateHelper::GetSimplePEPredicateExtractor(const conj_graph_pack& gp,
                                                                         const std::string path_to_contigs,
                                                                         const size_t length_threshold) {
    path_extend::validation::FilteredReferencePathHelper helper(gp);
    auto contig_paths = helper.GetFilteredReferencePathsFromLength(path_to_contigs,
                                                                   length_threshold);
    path_extend::validation::StrictTransitionStorageBuilder transition_builder;
    auto transition_storage = transition_builder.GetTransitionStorage(contig_paths);
    INFO("PE transition storage size: " << transition_storage.size());
    SimplePEPredicate pe_score_extractor(transition_storage);
    return pe_score_extractor;
}
}
