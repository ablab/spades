//
// Created by andrey on 21.09.15.
//

#ifndef PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
#define PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP

#include "scaffold_graph.hpp"
#include "visualization/graph_colorer.hpp"
#include "visualization/graph_labeler.hpp"
#include "visualization/graph_printer.hpp"

namespace path_extend {

namespace scaffold_graph {

using namespace visualization;


class ScaffoldGraphLabeler : public graph_labeler::GraphLabeler<ScaffoldGraph> {

private:
    const ScaffoldGraph &graph_;

    const std::map<VertexId, std::string> &additional_vertex_labels_;

public:
    ScaffoldGraphLabeler(const ScaffoldGraph &graph, const std::map<VertexId, std::string> &additional_vertex_labels):
        graph_(graph), additional_vertex_labels_(additional_vertex_labels) {
    }

    std::string label(VertexId v) const;

    std::string label(EdgeId e) const;
};


class ScaffoldEdgeColorer : public graph_colorer::ElementColorer<ScaffoldGraph::EdgeId> {
private:
    static const std::map<size_t, std::string> color_map;

    static const std::string default_color;

public:
    std::string GetValue(ScaffoldGraph::EdgeId e) const;
};


class ScaffoldVertexSetColorer : public graph_colorer::ElementColorer<ScaffoldGraph::VertexId> {
private:
    std::set<ScaffoldGraph::VertexId> vertex_set_;

public:
    ScaffoldVertexSetColorer(const std::set<ScaffoldGraph::VertexId> &vertex_set): vertex_set_(vertex_set) {
    }
    std::string GetValue(ScaffoldGraph::VertexId v) const;
};

class ScaffoldGraphVisualizer {
private:
    const ScaffoldGraph &graph_;

    const std::map<ScaffoldGraph::VertexId, std::string> &additional_vertex_labels_;

private:
    void Visualize(graph_printer::GraphPrinter<ScaffoldGraph> &printer);

public:
    ScaffoldGraphVisualizer(const ScaffoldGraph &graph,
                            const std::map<ScaffoldGraph::VertexId, std::string> &additional_vertex_labels) :
            graph_(graph),
            additional_vertex_labels_(additional_vertex_labels){
    }

    void Visualize(std::ostream &os, graph_colorer::CompositeGraphColorer<ScaffoldGraph> &colorer);
};

} //scaffold_graph
} //path_extend


#endif //PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
