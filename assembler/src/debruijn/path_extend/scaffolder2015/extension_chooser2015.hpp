//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "path_extend/extension_chooser.hpp"
#include "connection_condition2015.hpp"
#include "genome_consistance_checker.hpp"
#include "logger/logger.hpp"
#include <map>
#include <set>
namespace path_extend {
class ExtensionChooser2015: public ScaffoldingExtensionChooser {
private:
    shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges_;
// for possible connections e1 and e2 if weight(e1) > relative_weight_threshold_ * weight(e2) then e2 will be ignored
    double relative_weight_threshold_;
    PairedLibConnectionCondition paired_connection_condition_;
    AssemblyGraphConnectionCondition graph_connection_condition_;
// weight < absolute_weight_threshold_ will be ignored
    size_t absolute_weight_threshold_;
// multiplicator for the pairs which are connected in graph.
    double graph_connection_bonus_;

protected:
//If path contains no unique edges return -1
    pair<EdgeId, int> FindLastUniqueInPath(const BidirectionalPath& path) const;
//Find all possible next unique edges confirmed with mate-pair information. (absolute/relative)_weight_threshold_ used for filtering
    EdgeContainer FindNextUniqueEdge(const EdgeId from) const;
        DECL_LOGGER("ExtensionChooser2015")
public:
    ExtensionChooser2015(const Graph& g, shared_ptr<WeightCounter> wc, double is_scatter_coeff,
                         shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges ,double relative_threshold, size_t lib_index):
            ScaffoldingExtensionChooser(g, wc, is_scatter_coeff), unique_edges_(unique_edges), relative_weight_threshold_(relative_threshold), paired_connection_condition_(g,
            wc->get_libptr(), lib_index,
//TODO: constants are subject to reconsider
            0), graph_connection_condition_(g, 2*unique_edges_->GetMinLength()), absolute_weight_threshold_(2), graph_connection_bonus_(2) {
        INFO("ExtensionChooser2015 created");
    }
/* @param edges are really not used and left for compatibility
 * @returns possible next edge if there is unique one, else returns empty container
 *
 */

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override;
};


}
