#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_connection_conditions.hpp"

namespace path_extend {
    class SimplePEPredicate: public ScaffoldEdgePredicate {
        const path_extend::validation::ContigTransitionStorage transition_storage_;

     public:
        explicit SimplePEPredicate(const path_extend::validation::ContigTransitionStorage& transition_storage):
            transition_storage_(transition_storage) {}

        bool Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const override;
    };


    class SimplePEPredicateHelper {
     public:
        SimplePEPredicate GetSimplePEPredicateExtractor(const conj_graph_pack& gp, const std::string path_to_contigs,
                                                        const size_t length_threshold);
    };
}