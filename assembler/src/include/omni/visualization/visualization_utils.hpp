#pragma once

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


#include "graph_printer.hpp"
#include "omni/omni_utils.hpp"
#include "omni/dijkstra_tools/dijkstra_helper.hpp"
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
		shared_ptr<GraphSplitter<Graph>> inner_splitter,
		shared_ptr<GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph> &labeler) {
	EmptyGraphLinker<Graph> linker;
//	shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
    auto filter = make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CollectingSplitterWrapper<Graph>>(inner_splitter, filter);
	omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(*splitter, folder_name);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph>& gc,
		const string& file_name, shared_ptr<GraphColorer<Graph>> colorer,
		const GraphLabeler<Graph> &labeler) {
    EmptyGraphLinker<Graph> linker;
    BorderDecorator<Graph> component_colorer(gc, *colorer, "yellow");
    ofstream os;
    os.open(file_name);
	omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), true).Visualize(gc, os, labeler, component_colorer, linker);
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
    edge_colorer->AddColorer(make_shared<SetColorer<Graph>>(g, path.sequence(), "green"));
    shared_ptr<GraphColorer<Graph>> resulting_colorer =  make_shared<CompositeGraphColorer<Graph>>(colorer, edge_colorer);
    shared_ptr<GraphSplitter<Graph>> rs = ReliableSplitterAlongPath<Graph>(g, path);
    auto filter = make_shared<omnigraph::SmallComponentFilter<Graph>>(g, 3);
    shared_ptr<GraphSplitter<Graph>> splitter = make_shared<omnigraph::CondensingSplitterWrapper<Graph>>(rs, filter);
    WriteComponents<Graph>(g, folder_name, splitter, resulting_colorer, labeler);
}

//static void WriteFilteredComponents(const Graph& g,
//		const string& folder_name,
//		shared_ptr<GraphComponentFilter<Graph>> filter,
//		shared_ptr<GraphSplitter<Graph>> splitter,
//		shared_ptr<GraphColorer<Graph>> colorer,
//		const GraphLabeler<Graph> &labeler) {
//	EmptyGraphLinker<Graph> linker;
////	shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
//	omnigraph::FilteringSplitterWrapper<Graph> filtered_splitter(splitter, filter);
//	omnigraph::visualization::SplittingGraphVisualizer<Graph>(g, labeler, *colorer, linker).SplitAndVisualize(filtered_splitter, folder_name);
//}

}
}
