//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "extension_chooser.hpp"
#include "connection_condition2015.hpp"
#include "genome_consistance_checker.hpp"
#include "logger/logger.hpp"
#include <map>
#include <set>
namespace path_extend {
class ExtensionChooser2015: public ScaffoldingExtensionChooser {
private:
    shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges_;
    double relative_weight_threshold_;
    PairedLibConnectionCondition paired_connection_condition_;
    AssemblyGraphConnectionCondition graph_connection_condition_;
// < absolute threshold will be ignored
    size_t absolute_weight_threshold_;
    double graph_connection_bonus_;

protected:
    pair<EdgeId, int> FindLastUniqueInPath(const BidirectionalPath& path) const;
    EdgeContainer FindNextUniqueEdge(const EdgeId from) const;
        DECL_LOGGER("ExtensionChooser2015")
public:
    ExtensionChooser2015(const Graph& g, shared_ptr<WeightCounter> wc, double priority, double is_scatter_coeff,
                         shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges ,double relative_threshold, size_t lib_index):
            ScaffoldingExtensionChooser(g, wc, is_scatter_coeff), unique_edges_(unique_edges), relative_weight_threshold_(relative_threshold), paired_connection_condition_(g,
            wc->get_libptr(), lib_index,
//TODO: constants
            0), graph_connection_condition_(g, 2*unique_edges_->GetMinLength()), absolute_weight_threshold_(2), graph_connection_bonus_(2) {
        INFO("ExtensionChooser2015 created");
    }
//edges are really not used and left for compatibility
    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override;
};


}
