#pragma once

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


#include "graph_printer.hpp"
#include "omni/omni_utils.hpp"
#include "omni/dijkstra.hpp"
#include "omni/splitters.hpp"
#include "omni/abstract_conjugate_graph.hpp"
#include "omni/abstract_nonconjugate_graph.hpp"
#include "omni/graph_component.hpp"
#include "visualizers.hpp"
#include "vertex_linker.hpp"

namespace omnigraph {
namespace visualization {


template<class Graph>
void WriteComponents(const Graph& g,
		const string& folder_name,
		shared_ptr<GraphSplitter<Graph>> splitter,
		shared_ptr<GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph> &labeler) {
	EmptyGraphLinker<Graph> linker;
//	shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
	omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(*splitter, folder_name);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph>& gc,
		const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    ofstream os;
    os.open(file_name);
	omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), true).Visualize(gc, os, labeler, *colorer, linker);
	os.close();
}

template<class Graph>
void WriteSimpleComponent(const GraphComponent<Graph>& gc,
		const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    ofstream os;
    os.open(file_name);
	omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), false).Visualize(gc, os, labeler, *colorer, linker);
	os.close();
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g, Path<typename Graph::EdgeId> path,
		const string& folder_name, shared_ptr<GraphColorer<Graph>> colorer, const GraphLabeler<Graph> &labeler) {
    auto edge_colorer = make_shared<CompositeEdgeColorer<Graph>>("black");
    edge_colorer->AddColorer(colorer);
    edge_colorer->AddColorer(make_shared<PathColorer<Graph>>(g, path, "green"));
    shared_ptr<GraphColorer<Graph>> resulting_colorer =  make_shared<CompositeGraphColorer<Graph>>(colorer, edge_colorer);
    shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitterAlongPath<Graph>(g, path);
    WriteComponents<Graph>(g, folder_name, splitter, resulting_colorer, labeler);
}
}
}
