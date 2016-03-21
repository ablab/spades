//
// Created by andrey on 21.09.15.
//

#ifndef PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
#define PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP

#include "pipeline/graphio.hpp"
#include "scaffold_graph.hpp"

namespace path_extend { namespace scaffold_graph {

using namespace omnigraph::visualization;


class ScaffoldGraphLabeler : public GraphLabeler<ScaffoldGraph> {

private:
    const ScaffoldGraph &graph_;

public:
    ScaffoldGraphLabeler(const ScaffoldGraph &graph) : graph_(graph) {
    }

    string label(VertexId v) const;

    string label(EdgeId e) const;
};


class ScaffoldEdgeColorer : public ElementColorer<ScaffoldGraph::EdgeId> {
private:
    static const map<size_t, string> color_map;

    static const string default_color;

public:
    string GetValue(ScaffoldGraph::EdgeId e) const;
};


class ScaffoldVertexSetColorer : public ElementColorer<ScaffoldGraph::VertexId> {
 private:
  set<ScaffoldGraph::VertexId> vertex_set_;

 public:
  ScaffoldVertexSetColorer(const set<ScaffoldGraph::VertexId>& vertex_set): vertex_set_(vertex_set) {
  }

    string GetValue(ScaffoldGraph::VertexId v) const;
};

class ScaffoldGraphVisualizer {

    const ScaffoldGraph &graph_;
    const bool paired_;

private:
    void Visualize(GraphPrinter<ScaffoldGraph> &printer);

public:
    ScaffoldGraphVisualizer(const ScaffoldGraph &graph, bool paired = true) :
            graph_(graph), paired_(paired) {
    }

    void Visualize(ostream &os, CompositeGraphColorer<ScaffoldGraph>& colorer);
};

} //scaffold_graph
} //path_extend


#endif //PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
