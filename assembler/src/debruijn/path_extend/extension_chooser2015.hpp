//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "extension_chooser.hpp"
#include "genome_consistance_checker.hpp"
#include "logger/logger.hpp"
#include <map>
#include <set>
namespace path_extend {
class ExtensionChooser2015: public ScaffoldingExtensionChooser {
private:
    shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges_;
    double relative_weight_threshold_;
protected:
    int CountMedian(std::vector<std::pair<int, double> >& histogram) const;
    void FindBestFittedEdges(const BidirectionalPath& path, const std::set<EdgeId>& edges, EdgeContainer& result) const;
    std::set<EdgeId> FindCandidates(const BidirectionalPath& path) const ;
    DECL_LOGGER("ExtensionChooser2015")
public:
    ExtensionChooser2015(const Graph& g, shared_ptr<WeightCounter> wc, double priority, double is_scatter_coeff,
                         shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges ,double relative_threshold):
            ScaffoldingExtensionChooser(g, wc, is_scatter_coeff), unique_edges_(unique_edges), relative_weight_threshold_(relative_threshold) {
        INFO("ExtensionChooser2015 created");
    }

    EdgeContainer Filter(const BidirectionalPath& path, const EdgeContainer& edges) const override {
        set<EdgeId> candidates = FindCandidates(path);
        EdgeContainer result;
        FindBestFittedEdges(path, candidates, result);
        return result;
    }
};


}
