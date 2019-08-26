
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by andrey on 21.09.15.
//

#ifndef PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
#define PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "visualization/graph_colorer.hpp"
#include "visualization/graph_printer.hpp"

namespace path_extend {

namespace scaffolder {

using namespace visualization;


  class ScaffoldGraphLabeler : public graph_labeler::GraphLabeler<scaffold_graph::ScaffoldGraph> {

private:
    const scaffold_graph::ScaffoldGraph &graph_;

    const std::map<VertexId, std::string> &additional_vertex_labels_;

public:
    ScaffoldGraphLabeler(const scaffold_graph::ScaffoldGraph &graph,
                         const std::map<VertexId, std::string> &additional_vertex_labels):
        graph_(graph), additional_vertex_labels_(additional_vertex_labels) {
    }

    std::string label(VertexId v) const;

    std::string label(EdgeId e) const;
};


class ScaffoldEdgeColorer : public graph_colorer::ElementColorer<scaffold_graph::ScaffoldGraph::EdgeId> {
private:
    static const std::map<size_t, std::string> color_map;

    static const std::string default_color;

public:
    std::string GetValue(scaffold_graph::ScaffoldGraph::EdgeId e) const;
};


class ScaffoldVertexSetColorer : public graph_colorer::ElementColorer<scaffold_graph::ScaffoldGraph::VertexId> {
private:
    std::set<scaffold_graph::ScaffoldGraph::VertexId> vertex_set_;

public:
    ScaffoldVertexSetColorer(const std::set<scaffold_graph::ScaffoldGraph::VertexId> &vertex_set): vertex_set_(vertex_set) {
    }
    std::string GetValue(scaffold_graph::ScaffoldGraph::VertexId v) const;
};

class ScaffoldGraphVisualizer {
private:
    const scaffold_graph::ScaffoldGraph &graph_;

    const std::map<scaffold_graph::ScaffoldGraph::VertexId, std::string> &additional_vertex_labels_;

private:
    void Visualize(graph_printer::GraphPrinter<scaffold_graph::ScaffoldGraph> &printer);

public:
    ScaffoldGraphVisualizer(const scaffold_graph::ScaffoldGraph &graph,
                            const std::map<scaffold_graph::ScaffoldGraph::VertexId, std::string> &additional_vertex_labels) :
            graph_(graph),
            additional_vertex_labels_(additional_vertex_labels){
    }

    void Visualize(std::ostream &os, graph_colorer::CompositeGraphColorer<scaffold_graph::ScaffoldGraph> &colorer);
};

} //scaffold_graph
} //path_extend


#endif //PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
