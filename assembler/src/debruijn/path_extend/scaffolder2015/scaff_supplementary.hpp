#pragma once
#include "graph_pack.hpp"
#include "logger/logger.hpp"

namespace path_extend {
    typedef debruijn_graph::EdgeId EdgeId;

/* Storage of presumably unique, relatively long edges. Filled by ScaffoldingUniqueEdgeAnalyzer
 *
 */
    class ScaffoldingUniqueEdgeStorage {
        friend class ScaffoldingUniqueEdgeAnalyzer;
    private:
        set <EdgeId> unique_edges_;
        size_t min_unique_length_;
    public:
        ScaffoldingUniqueEdgeStorage(): unique_edges_(){
            DEBUG("storage created, empty");
        }

        bool IsUnique(EdgeId e) const {
            return (unique_edges_.find(e) != unique_edges_.end());
        }

        decltype(unique_edges_.begin()) begin() const {
            return unique_edges_.begin();
        }

        decltype(unique_edges_.end()) end() const {
            return unique_edges_.end();
        }

        size_t size() const {
            return unique_edges_.size();
        }
        size_t GetMinLength() const {
            return min_unique_length_;
        }
        void SetMinLength(size_t min_length)  {
            min_unique_length_ = min_length;
        }

        const set<EdgeId>& GetSet() const {
            return unique_edges_;
        }
   protected:
        DECL_LOGGER("ScaffoldingUniqueEdgeStorage")

    };

/* Auxillary class required to fillin the unique edge storage.
 *
 */
    class ScaffoldingUniqueEdgeAnalyzer {

    ;
    private:
        const debruijn_graph::conj_graph_pack &gp_;
        size_t length_cutoff_;
        double median_coverage_;
        double relative_coverage_variation_;
    protected:
        DECL_LOGGER("ScaffoldingUniqueEdgeAnalyzer")


        void SetCoverageBasedCutoff();
    public:
        ScaffoldingUniqueEdgeAnalyzer(const debruijn_graph::conj_graph_pack &gp, size_t apriori_length_cutoff, double max_relative_coverage):gp_(gp), length_cutoff_(apriori_length_cutoff), relative_coverage_variation_(max_relative_coverage){
            SetCoverageBasedCutoff();
        }
        void FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage_);
    };
}


