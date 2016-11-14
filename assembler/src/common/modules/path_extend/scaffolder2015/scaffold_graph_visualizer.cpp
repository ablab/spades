//
// Created by andrey on 21.09.15.
//

#include "scaffold_graph_visualizer.hpp"

namespace path_extend{ namespace scaffold_graph {

const map<size_t, string> ScaffoldEdgeColorer::color_map =
        {{(size_t) -1, "black"},
         {0, "red"},
         {1, "blue"},
         {2, "green"},
         {3, "magenta"},
         {4, "orange"},
         {5, "cyan"}};

const string ScaffoldEdgeColorer::default_color = "black";

string ScaffoldGraphLabeler::label(EdgeId e) const {
    return "ID: " + ToString(e.getId()) +
        "\\n Weight: " + ToString(e.getWeight()) +
        "\\n Lib#: " + ToString(e.getColor());
}

string ScaffoldGraphLabeler::label(VertexId v) const {
    auto it = additional_vertex_labels_.find(v);
    string additional_label = it == additional_vertex_labels_.end() ? "" : it->second + "\n";
    return "ID: " + ToString(graph_.int_id(v)) +
        "\\n Len: " + ToString(graph_.AssemblyGraph().length(v)) +
        "\\n Cov: " + ToString(graph_.AssemblyGraph().coverage(v)) + "\n" +
        additional_label;
}

void ScaffoldGraphVisualizer::Visualize(graph_printer::GraphPrinter<ScaffoldGraph> &printer) {
    printer.open();
    printer.AddVertices(graph_.vbegin(), graph_.vend());
    for (const auto& e : graph_.edges()) {
        printer.AddEdge(e);
    }
    printer.close();
}

void ScaffoldGraphVisualizer::Visualize(ostream &os, graph_colorer::CompositeGraphColorer<ScaffoldGraph>& colorer) {
    ScaffoldGraphLabeler labeler(graph_, additional_vertex_labels_);
    vertex_linker::EmptyGraphLinker<ScaffoldGraph> linker;

    graph_printer::SingleGraphPrinter <ScaffoldGraph> printer(graph_, os, labeler, colorer, linker);
    Visualize(printer);
}

string ScaffoldEdgeColorer::GetValue(ScaffoldGraph::EdgeId e) const {
    auto it = color_map.find(e.getColor());
    if (it != color_map.end()) {
        return it->second;
    }
    return default_color;
}

string ScaffoldVertexSetColorer::GetValue(ScaffoldGraph::VertexId v) const {
    if (vertex_set_.count(v) > 0)
        return "white";
    return "yellow";
}
} //scaffold_graph
} //path_extend



