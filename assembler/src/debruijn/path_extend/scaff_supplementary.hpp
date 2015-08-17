//
// Created by lab42 on 8/18/15.
//

#ifndef PROJECT_PE_SCAFF_SUPPLEMENTARY_HPP
#define PROJECT_PE_SCAFF_SUPPLEMENTARY_HPP

#include "graph_pack.hpp"

namespace path_extend {
    typedef debruijn_graph::EdgeId EdgeId;
    class ScaffoldingGenomeRepresentation {


    };

    class ScaffoldingUniqueEdgeStorage {
        friend class ScaffoldingUniqueEdgeAnalyzer;
    private:
        set <EdgeId> unique_edges_;
    public:

        bool IsUnique(EdgeId e) {
            return (unique_edges_.find(e) != unique_edges_.end());
        }

        set<EdgeId>::iterator begin() {
            return unique_edges_.begin();
        }

        set<EdgeId>::iterator end() {
            return unique_edges_.end();
        }
    };


    class ScaffoldingUniqueEdgeAnalyzer {

    ;
    private:
        const debruijn_graph::conj_graph_pack &gp_;
        size_t length_cutoff_;
        double median_coverage_;
        double max_relative_coverage_;
    protected:
        DECL_LOGGER("ScaffoldingUniqueEdgeAnalyzer")

        void SetInsertSizeBasedCutoff() {
//TODO: here should be something interesting, like max(cutoff, max(IS))
            return;
        }

        void SetCoverageBasedCutoff();
    public:
//default: 1000, 1.5 ?
        ScaffoldingUniqueEdgeAnalyzer(const debruijn_graph::conj_graph_pack &gp, size_t apriori_length_cutoff, double max_relative_coverage):gp_(gp), length_cutoff_(apriori_length_cutoff), max_relative_coverage_(max_relative_coverage){
            SetInsertSizeBasedCutoff();
            SetCoverageBasedCutoff();
        }
        void FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage_);
    };
}

#endif //PROJECT_PE_SCAFF_SUPPLEMENTARY_HPP
