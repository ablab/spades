//
// Created by andrey on 21.09.15.
//

#ifndef PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
#define PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP

#include "graphio.hpp"
#include "scaffold_graph.hpp"

namespace path_extend { namespace scaffold_graph {

using namespace omnigraph::visualization;


class ScaffoldGraphLabeler : public GraphLabeler<ScaffoldGraph> {

private:
    const ScaffoldGraph &graph_;

public:
    ScaffoldGraphLabeler(const ScaffoldGraph &graph) : graph_(graph) {

    }

    string label(VertexId v) const {
        return "ID: " + ToString(graph_.int_id(v)) +
                        "\\n Len: " + ToString(graph_.AssemblyGraph().length(v)) +
                        "\\n Cov: " + ToString(graph_.AssemblyGraph().coverage(v));
    }

    string label(EdgeId e) const {
        return "ID: " + ToString(e.getId()) +
                        "\\n Weight: " + ToString(e.getWeight()) +
                        "\\n Lib#: " + ToString(e.getColor());
    }
};

class ScaffoldEdgeColorer : public ElementColorer<ScaffoldGraph::EdgeId> {
private:
    static const map<size_t, string> color_map;

    static const string default_color;

public:
    string GetValue(ScaffoldGraph::EdgeId e) const {
        auto it = color_map.find(e.getColor());
        if (it != color_map.end()) {
            return it->second;
        }
        return default_color;
    }
};

class ScaffoldVertexSetColorer : public ElementColorer<ScaffoldGraph::VertexId> {
 private:
  set<ScaffoldGraph::VertexId> vertex_set_;

 public:
  ScaffoldVertexSetColorer(const set<ScaffoldGraph::VertexId>& vertex_set): vertex_set_(vertex_set) {
  }

  string GetValue(ScaffoldGraph::VertexId v) const {
      if (vertex_set_.count(v) > 0)
          return "white";
      return "yellow";
  }
};



class ScaffoldGraphVisualizer {

    const ScaffoldGraph &graph_;
    const bool paired_;

private:
    void Visualize(GraphPrinter<ScaffoldGraph> &printer) {
        printer.open();
        printer.AddVertices(graph_.vbegin(), graph_.vend());
        for (auto e_it = graph_.ebegin(); e_it != graph_.eend(); ++e_it) {
            printer.AddEdge(*e_it);
        }
        printer.close();
    }

public:
    ScaffoldGraphVisualizer(const ScaffoldGraph &graph, bool paired = true) :
            graph_(graph), paired_(paired) {
    }

    void Visualize(ostream &os, CompositeGraphColorer<ScaffoldGraph>& colorer) {
        ScaffoldGraphLabeler labeler(graph_);
        EmptyGraphLinker<ScaffoldGraph> linker;

        if (paired_) {
            PairedGraphPrinter <ScaffoldGraph> printer(graph_, os, labeler, colorer, linker);
            Visualize(printer);
        } else {
            SingleGraphPrinter <ScaffoldGraph> printer(graph_, os, labeler, colorer, linker);
            Visualize(printer);
        }
    }
};


} //scaffold_graph
} //path_extend


#endif //PROJECT_SCAFFOLD_GRAPH_VISUALIZER_HPP
