#pragma once

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


#include "standard_base.hpp"
#include "graph_printer.hpp"
namespace omnigraph {
namespace visualization {

//DECL_LOGGER("omg.gvis")

template<class Graph>
class ComponentVisualizer {
	const Graph& graph_;
	const bool paired_;

private:
	void Visualize(const GraphComponent<Graph>& component, GraphPrinter<Graph> &printer) {
		printer.open();
		printer.AddVertices(component.vertices().begin(), component.vertices().end());
		for (auto e_it = component.e_begin(); e_it != component.e_end();
				++e_it) {
			printer.AddEdge(*e_it);
		}
		printer.close();
	}

public:
	ComponentVisualizer(const Graph& graph, bool paired = true) :
		graph_(graph), paired_(paired) {
	}

	void Visualize(const GraphComponent<Graph>& component, ostream &os,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker) {
		if(paired_) {
			PairedGraphPrinter<Graph> printer(graph_, os, labeler, colorer, linker);
			Visualize(component, printer);
		} else {
			SingleGraphPrinter<Graph> printer(graph_, os, labeler, colorer, linker);
			Visualize(component, printer);
		}
	}

	void Visualize(ostream &os,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker) {
		GraphComponent<Graph> component(graph_, graph_.begin(), graph_.end(), false);
		Visualize(component, os, labeler, colorer, linker);
	}
};


template<class Graph>
class ComponentNameGenerator {
public:
	virtual string ComponentName(const GraphComponent<Graph>& component) = 0;

	virtual ~ComponentNameGenerator() {
	}
};

template<class Graph>
class SimpleCountingComponentNameGenerator: public ComponentNameGenerator<Graph> {
private:
	string name_;
	string extension_;
	size_t cnt_;
public:
	SimpleCountingComponentNameGenerator(string name, string extension): name_(name), extension_(extension), cnt_(0) {
	}

	string ComponentName(const GraphComponent<Graph>& component) {
		cnt_++;
		stringstream ss;
		ss << name_ << "_" << cnt_;
		if(component.name().size() > 0)
			ss << "_" << component.name();
		ss << "." << extension_;
		return ss.str();
	}
};

template<class Graph>
class SplittingGraphVisualizer {
private:
	static const size_t DEFAULT_MAX_COMPONENT_NUMBER = 500;
	const Graph& graph_;
	const GraphLabeler<Graph> &labeler_;
	const GraphColorer<Graph> &colorer_;
	const VertexLinker<Graph> &linker_;
	const bool paired_;
	const size_t max_component_number_;

    string ComponentFileName(size_t cnt, const string &folder, const GraphComponent<Graph>& component) {
        stringstream ss;
        ss << folder << cnt;
        if(component.name().size() > 0)
            ss << "graph_" << component.name();
        ss << ".dot";
        return ss.str();
    }

public:
	SplittingGraphVisualizer(const Graph& graph,
			const GraphLabeler<Graph> &labeler,
			const GraphColorer<Graph> &colorer,
			const VertexLinker<Graph> &linker,
			bool paired = true,
			size_t max_component_number = DEFAULT_MAX_COMPONENT_NUMBER) :
			graph_(graph), labeler_(labeler), colorer_(colorer), linker_(linker), paired_(paired), max_component_number_(max_component_number) {
	}

	size_t SplitAndVisualize(GraphSplitter<Graph> &splitter, const string &folder) {
	    INFO("Writing components to folder " << folder);
		ComponentVisualizer<Graph> visualizer(graph_, paired_);
		size_t cnt = 0;
		while(splitter.HasNext()) {
			if(cnt > max_component_number_) {
				INFO("The number of graph components exceeded " << max_component_number_ << ". Aborting current visualization.");
				break;
			}
			cnt++;
			GraphComponent<Graph> component = splitter.Next();
			BorderDecorator<Graph> border_colorer(component, colorer_, "yellow");
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

