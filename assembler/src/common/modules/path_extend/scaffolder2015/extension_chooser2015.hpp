//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "modules/path_extend/extension_chooser.hpp"
#include "connection_condition2015.hpp"
#include "modules/genome_consistance_checker.hpp"
#include "utils/logger/logger.hpp"
#include <map>
#include <set>
namespace path_extend {

class ExtensionChooser2015: public ScaffoldingExtensionChooser {
    static const int MIN_N_QUANTITY = 10;
    shared_ptr<ConnectionCondition> lib_connection_condition_;
    const ScaffoldingUniqueEdgeStorage& unique_edges_;

    // for possible connections e1 and e2 if weight(e1) > relative_weight_threshold_ * weight(e2) then e2 will be ignored
    double relative_weight_threshold_;
    AssemblyGraphConnectionCondition graph_connection_condition_;
    // weight < absolute_weight_threshold_ will be ignored
    size_t absolute_weight_threshold_;
    // multiplicator for the pairs which are connected in graph.
    double graph_connection_bonus_;
    bool use_graph_connectivity_;

    //If path contains no unique edges return -1
    pair<EdgeId, int> FindLastUniqueInPath(const BidirectionalPath& path) const;
    //Find all possible next unique edges confirmed with mate-pair information. (absolute/relative)_weight_threshold_ used for filtering
    EdgeContainer FindNextUniqueEdge(const EdgeId from) const;
public:
    ExtensionChooser2015(const Graph& g,
                         shared_ptr<WeightCounter> wc,
                         shared_ptr<ConnectionCondition> condition,
                         const ScaffoldingUniqueEdgeStorage& unique_edges,
                         double cl_weight_threshold,
                         double is_scatter_coeff,
                         double relative_threshold,
                         bool use_graph_connectivity = true):
            //TODO: constants are subject to reconsider
            ScaffoldingExtensionChooser(g, wc, cl_weight_threshold, is_scatter_coeff),
            lib_connection_condition_(condition),
            unique_edges_(unique_edges),
            relative_weight_threshold_(relative_threshold),
            graph_connection_condition_(g, 2 * unique_edges_.min_length(), unique_edges),
            //TODO to config!
            absolute_weight_threshold_(2),
            graph_connection_bonus_(2),
            use_graph_connectivity_(use_graph_connectivity) {
        INFO("ExtensionChooser2015 created");
    }

    /* @param edges are really not used and left for compatibility
     * @returns possible next edge if there is unique one, else returns empty container
     */
    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer&) const override;
    void InsertAdditionalGaps(ExtensionChooser::EdgeContainer& result) const;

private:
    DECL_LOGGER("ExtensionChooser2015");
};


}
