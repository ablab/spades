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

    const map<VertexId, string>& additional_vertex_labels_;

public:
    ScaffoldGraphLabeler(const ScaffoldGraph &graph, const map<VertexId, string>& additional_vertex_labels):
        graph_(graph), additional_vertex_labels_(additional_vertex_labels) {
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
private:
    const ScaffoldGraph &graph_;

    const map<ScaffoldGraph::VertexId, string>& additional_vertex_labels_;

private:
    void Visualize(GraphPrinter<ScaffoldGraph> &printer);

public:
    ScaffoldGraphVisualizer(const ScaffoldGraph &graph,
                            const map<ScaffoldGraph::VertexId, string>& additional_vertex_labels) :
            graph_(graph),
            additional_vertex_labels_(additional_vertex_labels){
    }

    void Visualize(ostream &os, CompositeGraphColorer<ScaffoldGraph>& colorer);
};

} //scaffold_graph
} //path_extend


#endif //PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
