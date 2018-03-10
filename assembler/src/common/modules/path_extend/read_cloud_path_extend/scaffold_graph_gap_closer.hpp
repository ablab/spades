#pragma once

#include "scaffolder2015/scaffold_vertex.hpp"
#include "scaffolder2015/scaffold_graph.hpp"
namespace path_extend {

class TipFinder {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

 public:
    virtual ScaffoldVertex FindTip(const ScaffoldVertex& vertex) = 0;
};

class ScaffoldGraphGapCloser {
 protected:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
 public:
    virtual void CloseGaps(ScaffoldGraph& graph) const = 0;
};

class TipFinderGapCloser: public ScaffoldGraphGapCloser {
 protected:
    using ScaffoldGraphGapCloser::ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    shared_ptr<TipFinder> tip_finder_;

 public:
    void CloseGaps(ScaffoldGraph &graph) const override;
};
}