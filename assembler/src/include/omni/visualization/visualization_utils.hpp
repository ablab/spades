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
static set<typename Graph::VertexId> EndVertices(const Graph & graph,
                                 const set<typename Graph::EdgeId> &error_edges) {
    set<typename Graph::VertexId> result;
    for (auto it = error_edges.begin(); it != error_edges.end(); ++it) {
        result.insert(graph.EdgeEnd(*it));
        result.insert(graph.EdgeStart(*it));
    }
    return result;
}

template<class Graph>
void WriteErrors(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& folder_name,
		const Path<typename Graph::EdgeId> &path1,
		const Path<typename Graph::EdgeId> &path2) {
	shared_ptr<ElementColorer<typename Graph::EdgeId>> colorer =
            make_shared<CompositeEdgeColorer<Graph>>(
                    make_shared<PathColorer<Graph>>(g, path1, "red"),
                    make_shared<PathColorer<Graph>>(g, path2, "blue"), "black");
	GraphComponent<Graph> all(g, g.begin(), g.end());
	set<typename Graph::VertexId> black = EndVertices(g, colorer->ColoredWith(all.edges().begin(), all.edges().end(), "black"));
	shared_ptr<GraphSplitter<Graph>> splitter = StandardSplitter(g, black);
	shared_ptr<GraphComponentFilter<Graph>> checker = make_shared<ComponentSizeFilter<Graph>>(g, 1500, 2, 300);
	WriteComponents(g, splitter, checker, folder_name, *DefaultColorer(g, colorer), labeler);
}

//todo strange similar looking methods!!!
template<class Graph>
void WriteComponents(const Graph& graph,
		GraphSplitter<Graph> &splitter,
		const string& folder_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler) {
	EmptyGraphLinker<Graph> linker;
	omnigraph::visualization::SplittingGraphVisualizer<Graph>(graph, labeler, colorer, linker).SplitAndVisualize(splitter, folder_name);
}

template<class Graph>
void WriteComponents(const Graph& g,
		shared_ptr<GraphSplitter<Graph>> inner_splitter,
		shared_ptr<GraphComponentFilter<Graph>> &checker,
		const string& folder_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler) {
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	WriteComponents<Graph>(g, splitter, folder_name, colorer, labeler);
//    WriteComponents<Graph>(g, inner_splitter, folder_name, colorer, labeler);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph>& gc,
		const string& file_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler) {
    EmptyGraphLinker<Graph> linker;
    ofstream os;
    os.open(file_name);
	omnigraph::visualization::ComponentVisualizer<Graph>(gc.g(), true).Visualize(gc, os, labeler, colorer, linker);
	os.close();
}

template<class Graph>
void WriteComponents(const Graph& g, size_t split_edge_length,
		const string& folder_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler) {
	shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitter<Graph>(g, split_edge_length, 60);
	shared_ptr<GraphComponentFilter<Graph>> filter = make_shared<ComponentSizeFilter<Graph>>(g, split_edge_length, 2, 300);
	WriteComponents<Graph>(g, splitter, filter, folder_name, colorer, labeler);
}

//todo alert!!! magic constants!!!
//todo refactoring of params needed
template<class Graph>
void WriteComponentsAlongPath(const Graph& g, const AbstractFilter<vector<typename Graph::VertexId>> &checker,
		const GraphLabeler<Graph>& labeler, const string& file_name,
		size_t split_edge_length, size_t component_vertex_number,
		const MappingPath<typename Graph::EdgeId>& path,
		const GraphColorer<Graph>& colorer) {
	shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitterAlongPath<Graph>(g, path, split_edge_length);
	WriteComponents<Graph>(g, *splitter, checker, file_name,
			colorer, labeler);
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g,
                              const GraphLabeler<Graph>& labeler,
                              const string& folder_name, size_t split_edge_length,
                              size_t component_vertex_number,
                              Path<typename Graph::EdgeId> path,
                              const GraphColorer<Graph>& colorer) {
    shared_ptr<GraphSplitter<Graph>> splitter = ReliableSplitterAlongPath<Graph>(g, path, split_edge_length);
    WriteComponents<Graph>(g, *splitter, folder_name, colorer, labeler);
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g,
		const GraphLabeler<Graph>& labeler, const string& folder_name,
		size_t split_edge_length, size_t component_vertex_number,
		Path<typename Graph::EdgeId> path,
		Path<typename Graph::EdgeId> color1 = Path<typename Graph::EdgeId>(),
		Path<typename Graph::EdgeId> color2 = Path<typename Graph::EdgeId>()) {
    auto edge_colorer = make_shared<CompositeEdgeColorer<Graph>>("black");
    edge_colorer->AddColorer(make_shared<PathColorer<Graph>>(g, color1, "red"));
    edge_colorer->AddColorer(make_shared<PathColorer<Graph>>(g, color2, "blue"));
    edge_colorer->AddColorer(make_shared<PathColorer<Graph>>(g, path, "green"));
	shared_ptr<GraphColorer<Graph>> colorer =  DefaultColorer<Graph>(g, edge_colorer);
	WriteComponentsAlongPath(g, labeler, folder_name, split_edge_length, component_vertex_number, path, *colorer);
}

//todo alert!!! magic constants!!!
template<class Graph>
void WriteComponentsAlongGenome(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& folder_name,
		size_t split_edge_length,
		MappingPath<typename Graph::EdgeId> color1 = MappingPath<
				typename Graph::EdgeId>(),
		MappingPath<typename Graph::EdgeId> color2 = MappingPath<
				typename Graph::EdgeId>()) {
	WriteComponentsAlongPath<Graph>(g, labeler, folder_name, split_edge_length, 60,
			color1.simple_path(), color1.simple_path(), color2.simple_path());
}

template<class Graph>
void WriteComponentsAroundEdge(const Graph& g, typename Graph::EdgeId e,
		const string& file_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler, size_t length_bound = 40000, size_t max_size = 30) {
    shared_ptr<GraphSplitter<Graph>> splitter = EdgeNeighborhoodFinder<Graph>(g, e, max_size, length_bound);
	WriteComponents(g, *splitter, file_name, colorer, labeler);
}

template<class Graph>
void WriteComponentsAroundEdge(const Graph& g, const AbstractFilter<vector<typename Graph::VertexId>> &filter, typename Graph::EdgeId e,
		const string& file_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler) {
	shared_ptr<GraphSplitter<Graph>> splitter = EdgeNeighborhoodFinder<Graph>(g, e, 30, 40000);
	FilteringSplitterWrapper<Graph> filtering_splitter(splitter, filter);
	WriteComponents(g, filtering_splitter/*, "locality_of_edge_" + ToString(g_.int_id(edge))*/
	, file_name, colorer, labeler);
}

template<class Graph>
void WritePaired(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const Path<typename Graph::EdgeId> &path1 = Path<typename Graph::EdgeId>(),
		const Path<typename Graph::EdgeId> &path2 = Path<typename Graph::EdgeId>()) {
    WriteWholeGraph(g, true, labeler, file_name, path1, path2);
}

template<class Graph>
void WriteSimple(const Graph& g, const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const Path<typename Graph::EdgeId> &path1 = Path<typename Graph::EdgeId>(),
		const Path<typename Graph::EdgeId> &path2 = Path<typename Graph::EdgeId>()) {
    WriteWholeGraph(g, false, labeler, file_name, path1, path2);
}

template<class Graph>
void WriteWholeGraph(
        const Graph& g,
        bool paired,
        const GraphLabeler<Graph>& labeler,
        const string& file_name,
        const Path<typename Graph::EdgeId> &path1 = Path<typename Graph::EdgeId>(),
        const Path<typename Graph::EdgeId> &path2 = Path<typename Graph::EdgeId>()) {
    typedef typename Graph::EdgeId EdgeId;
    ofstream filestr(file_name);
    auto colorer = DefaultColorer(g, path1, path2);
    EmptyGraphLinker<Graph> linker;
    visualization::ComponentVisualizer<Graph>(g ,paired).Visualize(filestr, labeler, *colorer, linker);
    filestr.close();
}
}
}
