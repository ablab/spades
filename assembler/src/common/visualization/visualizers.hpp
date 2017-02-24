#pragma once

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************


#include "utils/standard_base.hpp"
#include "graph_printer.hpp"
#include <fstream>

using namespace omnigraph;

namespace visualization {

namespace visualizers {

//DECL_LOGGER("omg.gvis")

template<class Graph>
class ComponentVisualizer {
    const Graph &graph_;
    const bool paired_;

private:
    void Visualize(const GraphComponent<Graph> &component, graph_printer::GraphPrinter <Graph> &printer) {
        printer.open();
        printer.AddVertices(component.vertices().begin(), component.vertices().end());
        for (auto e_it = component.e_begin(); e_it != component.e_end();
             ++e_it) {
            printer.AddEdge(*e_it);
        }
        printer.close();
    }

public:
    ComponentVisualizer(const Graph &graph, bool paired = true) :
            graph_(graph), paired_(paired) {
    }

    void Visualize(const GraphComponent<Graph> &component, ostream &os,
                   const graph_labeler::GraphLabeler<Graph> &labeler,
                   const graph_colorer::GraphColorer<Graph> &colorer,
                   const vertex_linker::VertexLinker<Graph> &linker) {
        if (paired_) {
            graph_printer::PairedGraphPrinter<Graph> printer(graph_, os, labeler, colorer, linker);
            Visualize(component, printer);
        } else {
            graph_printer::SingleGraphPrinter<Graph> printer(graph_, os, labeler, colorer, linker);
            Visualize(component, printer);
        }
    }

    void Visualize(ostream &os,
                   const graph_labeler::GraphLabeler<Graph> &labeler,
                   const graph_colorer::GraphColorer<Graph> &colorer,
                   const vertex_linker::VertexLinker<Graph> &linker) {
        Visualize(GraphComponent<Graph>::WholeGraph(graph_), os, labeler, colorer, linker);
    }
};


template<class Graph>
class ComponentNameGenerator {
public:
    virtual string ComponentName(const GraphComponent<Graph> &component) = 0;

    virtual ~ComponentNameGenerator() {
    }
};

template<class Graph>
class SimpleCountingComponentNameGenerator : public ComponentNameGenerator<Graph> {
private:
    string name_;
    string extension_;
    size_t cnt_;
public:
    SimpleCountingComponentNameGenerator(string name, string extension) : name_(name), extension_(extension),
                                                                          cnt_(0) {
    }

    string ComponentName(const GraphComponent<Graph> &component) {
        cnt_++;
        stringstream ss;
        ss << name_ << "_" << cnt_;
        if (component.name().size() > 0)
            ss << "_" << component.name();
        ss << "." << extension_;
        return ss.str();
    }
};

template<class Graph>
class CountingSizeComponentNameGenerator : public ComponentNameGenerator<Graph> {
private:
    string name_;
    string extension_;
    size_t cnt_;
public:
    CountingSizeComponentNameGenerator(string name, string extension) : name_(name), extension_(extension),
                                                                        cnt_(0) {
    }

    string ComponentName(const GraphComponent<Graph> &component) {
        cnt_++;
        stringstream ss;
        ss << name_ << "_" << cnt_;
        if (component.name().size() > 0)
            ss << "_" << component.name();
        ss << "_size_" << component.size();
        ss << "." << extension_;

        return ss.str();
    }
};


template<class Graph>
class SplittingGraphVisualizer {
private:
    const Graph &graph_;
    const graph_labeler::GraphLabeler <Graph> &labeler_;
    const graph_colorer::GraphColorer <Graph> &colorer_;
    const vertex_linker::VertexLinker <Graph> &linker_;
    const bool paired_;
    const size_t max_component_number_;
    static const size_t DEFAULT_MAX_COMPONENT_NUMBER = 500;

    string ComponentFileName(size_t cnt, const string &folder, const GraphComponent<Graph> &component) {
        stringstream ss;
        ss << folder << cnt;
        if (component.name().size() > 0)
            ss << "graph_" << component.name();
        ss << ".dot";
        return ss.str();
    }

public:
    SplittingGraphVisualizer(const Graph &graph,
                             const graph_labeler::GraphLabeler <Graph> &labeler,
                             const graph_colorer::GraphColorer <Graph> &colorer,
                             const vertex_linker::VertexLinker <Graph> &linker,
                             bool paired = true,
                             size_t max_component_number = DEFAULT_MAX_COMPONENT_NUMBER) :
            graph_(graph), labeler_(labeler), colorer_(colorer), linker_(linker), paired_(paired),
            max_component_number_(max_component_number) {
    }

    size_t SplitAndVisualize(GraphSplitter<Graph> &splitter, const string &folder) {
        INFO("Writing components to folder " << folder);
        ComponentVisualizer<Graph> visualizer(graph_, paired_);
        size_t cnt = 0;
        while (splitter.HasNext()) {
            if (cnt > max_component_number_) {
                INFO("The number of graph components exceeded " << max_component_number_
                                                                << ". Aborting current visualization.");
                break;
            }
            cnt++;
            GraphComponent<Graph> component = splitter.Next();
            graph_colorer::BorderDecorator<Graph> border_colorer(component, colorer_, "yellow");
            ofstream os(ComponentFileName(cnt, folder, component));
            visualizer.Visualize(component, os, labeler_, border_colorer, linker_);
            os.close();
        }
        return cnt;
    }

private:
    DECL_LOGGER("SplittingGraphVisualizer");
};

}
}


