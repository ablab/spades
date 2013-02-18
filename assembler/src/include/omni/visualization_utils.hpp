//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "graph_printer.hpp"
#include "omni_utils.hpp"
#include "dijkstra.hpp"
#include "splitters.hpp"
#include "abstract_conjugate_graph.hpp"
#include "abstract_nonconjugate_graph.hpp"
#include "graph_component.hpp"

namespace omnigraph {

//DECL_LOGGER("omg.gvis")

using omnigraph::GraphLabeler;
using omnigraph::EmptyGraphLabeler;
using omnigraph::SmartEdgeIterator;
using omnigraph::Path;
using omnigraph::UnorientedDijkstra;

template<class Graph>
class GraphVisualizer {
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
protected:
	const Graph& g() {
		return g_;
	}
public:
	GraphVisualizer(const Graph& g) :
			g_(g) {

	}

	virtual ~GraphVisualizer() {

	}

	virtual void Visualize() = 0;

};

template<class Graph>
class PartialGraphVisualizer {
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;

protected:
	const Graph& g() {
		return g_;
	}

public:
	PartialGraphVisualizer(const Graph& g) :
			g_(g) {
	}

	virtual ~PartialGraphVisualizer() {

	}

	virtual void open() = 0;
	virtual void close() = 0;

	virtual void Visualize(const GraphComponent<Graph>& component) = 0;

};

template<class Graph>
class SimpleGraphVisualizer: public GraphVisualizer<Graph> {
	typedef GraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	GraphPrinter<Graph>& gp_;
public:
	SimpleGraphVisualizer(const Graph& g, GraphPrinter<Graph>& gp) :
			super(g), gp_(gp) {
	}

	virtual void Visualize() {
		gp_.open();
		TRACE("Visualize started");
		for (auto it = this->g().SmartVertexBegin(); !it.IsEnd(); ++it) {
			gp_.AddVertex(*it);
		}TRACE("Vertices printed");
		for (auto it = this->g().SmartEdgeBegin(); !it.IsEnd(); ++it) {
			gp_.AddEdge(*it);
		}TRACE("Edges printed");
		gp_.close();
	}

private:
	DECL_LOGGER("SimpleGraphVisualizer")
};

template<class Graph>
class AdapterGraphVisualizer: public GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef GraphVisualizer<Graph> super;
	PartialGraphVisualizer<Graph>& partial_visualizer_;
public:
	AdapterGraphVisualizer(const Graph& g,
			PartialGraphVisualizer<Graph>& partial_visualizer) :
			super(g), partial_visualizer_(partial_visualizer) {

	}

	virtual void Visualize() {
		partial_visualizer_.open();
		partial_visualizer_.Visualize(
				GraphComponent<Graph>(this->g()));
		partial_visualizer_.close();
	}
};

//todo rename!!!
template<class Graph>
class ColoredGraphVisualizer: public PartialGraphVisualizer<Graph> {
	typedef PartialGraphVisualizer<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	GraphPrinter<Graph>& printer_;
public:

	ColoredGraphVisualizer(const Graph& g, GraphPrinter<Graph>& printer) :
			base(g), printer_(printer) {
	}

	void open() {
		printer_.open();
	}

	void close() {
		printer_.close();
	}

	virtual void Visualize(const GraphComponent<Graph>& component) {
		printer_.AddVertices(component.vertices());
		for (auto e_it = component.e_begin(); e_it != component.e_end();
				++e_it) {
			printer_.AddEdge(*e_it);
		}
	}
};

//template<class Graph>
//void WriteToDotFile(const Graph& g, const GraphLabeler<Graph>& labeler,
//		const string& file_name, const string& graph_name = "my_graph") {
//	fstream filestr;
//	filestr.open(file_name.c_str(), fstream::out);
//	auto_ptr<GraphColorer<Graph>> colorer(DefaultColorer(g));
//	DotGraphPrinter<Graph> gpr(g, labeler, *colorer, graph_name, filestr);
//	SimpleGraphVisualizer<Graph> sgv(g, gpr);
//	sgv.Visualize();
//	filestr.close();
//}

//template<class Graph>
//void WriteSimple(const Graph& g, const GraphLabeler<Graph>& labeler,
//		const string& file_name, const string& graph_name) {
//	fstream filestr;
//	string simple_file_name(file_name);
//	//	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
//	filestr.open((simple_file_name).c_str(), fstream::out);
//	auto_ptr<GraphColorer<Graph>> colorer(DefaultColorer(g));
//	DotGraphPrinter<Graph> gp(g, labeler, *colorer, graph_name, filestr);
//	ColoredGraphVisualizer<Graph> gv(g, gp);
//	AdapterGraphVisualizer<Graph> result_vis(g, gv);
//	result_vis.Visualize();
//	filestr.close();
//}

template<class Graph>
void WriteSimple(const Graph& g, const GraphLabeler<Graph>& labeler,
		const string& file_name, const string& graph_name = "my_graph",
		const Path<typename Graph::EdgeId> &path1 = Path<typename Graph::EdgeId>(),
		const Path<typename Graph::EdgeId> &path2 = Path<typename Graph::EdgeId>()) {
	ofstream filestr(file_name);
	CompositeGraphColorer<Graph> colorer(new FixedColorer<typename Graph::VertexId>("white")
			, new MapColorer<typename Graph::EdgeId>(PathColorer<Graph>(g, path1, path2).ColorPath(), ""));
	DotGraphPrinter<Graph> gp(g, labeler, colorer, graph_name, filestr);
	SimpleGraphVisualizer<Graph> gv(g, gp);
	gv.Visualize();
	filestr.close();
}

template<class Graph>
void WritePaired(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const string& graph_name,
		const Path<typename Graph::EdgeId> &path1 = Path<typename Graph::EdgeId>(),
		const Path<typename Graph::EdgeId> &path2 = Path<typename Graph::EdgeId>()) {
	typedef typename Graph::EdgeId EdgeId;
	ofstream filestr(file_name);
	CompositeGraphColorer<Graph> colorer(new FixedColorer<typename Graph::VertexId>("")
			, new MapColorer<typename Graph::EdgeId>(PathColorer<Graph>(g, path1, path2).ColorPath(), ""));
	DotPairedGraphPrinter<Graph> gp(g, labeler, colorer, graph_name, filestr);
	SimpleGraphVisualizer<Graph> gv(g, gp);
	gv.Visualize();
	filestr.close();
}

template<class Graph>
class VisualizerFactory {
public:
	virtual auto_ptr<GraphPrinter<Graph>> GetPrinterInstance(
			const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer,
			const string &graph_name, ostream &os) = 0;

	virtual auto_ptr<PartialGraphVisualizer<Graph>> GetVisualizerInstance(GraphPrinter<Graph> &printer) = 0;

	virtual ~VisualizerFactory() {
	}
};

template<class Graph>
class ColoredVisualizerFactory: public VisualizerFactory<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const Graph &graph_;

	//todo dirty hack!!!
	auto_ptr<GraphPrinter<Graph>> PrinterInstance(
			const AbstractConjugateGraph<typename Graph::DataMaster>& hack,
			const Graph& graph,
			const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer,
			const string &graph_name,
			ostream &os) {
		return auto_ptr<GraphPrinter<Graph>>(new DotPairedGraphPrinter<Graph>(
				graph, labeler, colorer, graph_name, os));
	}

	auto_ptr<GraphPrinter<Graph>> PrinterInstance(
			const AbstractNonconjugateGraph<typename Graph::DataMaster>& hack,
			const Graph& graph,
			const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer,
			const string &graph_name,
			ostream &os) {
		return auto_ptr<GraphPrinter<Graph>>(new DotGraphPrinter<Graph>(graph, labeler, colorer, graph_name, os));
	}

public:
	ColoredVisualizerFactory(const Graph& graph) :
			graph_(graph) {
	}

	virtual auto_ptr<GraphPrinter<Graph>> GetPrinterInstance(
			const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer,
			const string &graph_name, ostream &os) {
		return PrinterInstance(graph_, graph_, labeler, colorer, graph_name, os);
	}

	virtual auto_ptr<PartialGraphVisualizer<Graph>> GetVisualizerInstance(GraphPrinter<Graph> &printer) {
		return auto_ptr<PartialGraphVisualizer<Graph>>(
				new ColoredGraphVisualizer<Graph>(graph_, printer));
	}

};

template<class Graph>
class ComponentGraphVisualizer: public GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	VisualizerFactory<Graph> &vis_factory_;
	ComponentSplitter<VertexId> &splitter_;
	const GraphLabeler<Graph>& labeler_;
	const GraphColorer<Graph>& colorer_;
	const string &file_name_;
	const string &graph_name_;
	size_t max_parts_number_;

	string ConstructComponentName(const string &file_name, size_t cnt,
			const string &component_name) {
		stringstream ss;
		ss << "_" << cnt << "_" << component_name;
		string res = file_name;
		res.insert(res.length() - 4, ss.str());
		return res;
	}

public:
	ComponentGraphVisualizer(const Graph &graph,
			VisualizerFactory<Graph> &factory,
			ComponentSplitter<typename Graph::VertexId> &splitter,
			const GraphLabeler<Graph>& labeler,
			const GraphColorer<Graph>& colorer,
			const string &file_name, const string &graph_name = "my_graph",
			size_t max_parts_number = 100) :
			GraphVisualizer<Graph>(graph), vis_factory_(factory), splitter_(
					splitter), labeler_(labeler), colorer_(colorer), file_name_(file_name), graph_name_(graph_name), max_parts_number_(
					max_parts_number) {
	}

	virtual ~ComponentGraphVisualizer() {
	}

	virtual void Visualize() {
		size_t cnt = 1;
		while (!splitter_.Finished() && cnt <= max_parts_number_) {
			auto vertices = splitter_.NextComponent();
			string component_name = ConstructComponentName(file_name_, cnt,
					splitter_.ComponentName());
			ofstream os;
			os.open(component_name.c_str());

			GraphComponent<Graph> component(this->g(), vertices.begin(), vertices.end());

			CompositeGraphColorer<Graph> local_colorer(
					new MapColorer<VertexId>(
							colorer_.GetColours(component.vertices())),
					new MapColorer<EdgeId>(
							colorer_.GetColours(component.edges())));
			auto_ptr<GraphPrinter<Graph>> printer(vis_factory_.GetPrinterInstance(labeler_,
							local_colorer, graph_name_, os));

			auto_ptr<PartialGraphVisualizer<Graph>> visualizer(vis_factory_.GetVisualizerInstance(*printer));

			visualizer->open();
			if (component.v_size() >= 1000) { // what the magic constant???
				WARN("Too large component " << component.v_size());
			}
			visualizer->Visualize(component);
			visualizer->close();
			os.close();
			cnt++;
		}
	}
};

template<class Graph>
string InsertComponentName(string file_name, string component) {
	string res = file_name;
	res.insert(res.length() - 4, "_" + component + "_");
	return res;
}

template<class Graph>
void WriteErrors(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		const string& graph_name,
		const Path<typename Graph::EdgeId> &path1/* = Path<typename Graph::EdgeId> ()*/,
		const Path<typename Graph::EdgeId> &path2/* = Path<typename Graph::EdgeId> ()*/) {
	PathColorer<Graph> path_colorer(g, path1, path2);
	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
	ErrorComponentSplitter<Graph> splitter(g, black);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	ComponentSizeFilter<Graph> checker(g, 1500, 2, 300);
	string error_file_name = InsertComponentName<Graph>(file_name, "error");
	WriteComponents(g, splitter, checker, error_file_name,
			*DefaultColorer(g, coloring), labeler);

//	PathColorer<Graph> path_colorer(g, path1, path2);
//	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
//	ErrorComponentSplitter<Graph> splitter(g, black);
//	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
//	ColoredVisualizerFactory<Graph> factory(g, labeler, coloring);
//	string error_file_name = InsertComponentName<Graph>(file_name, "error");
//	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, error_file_name,
//			graph_name, 200);
//	gv.Visualize();
	//	size_t cnt = 0;
	//	while (!splitter.Finished() && cnt < 100) {
	//		fstream filestr;
	//		filestr.open(ConstructComponentName<Graph> (file_name, cnt).c_str(),
	//				fstream::out);
	//		gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
	//		ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	//		auto component = splitter.NextComponent();
	//		gp.open();
	//		gv.Visualize(component);
	//		gp.close();
	//		cnt++;
	//	}
}

//todo strange similar looking methods!!!
template<class Graph>
void WriteComponents(const Graph& g,
		ComponentSplitter<typename Graph::VertexId> &splitter,
		const string& file_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler, const string& graph_name =
				"my_graph") {
	ColoredVisualizerFactory<Graph> factory(g);
	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, labeler, colorer, file_name,
			graph_name, 24000);
	gv.Visualize();
}

template<class Graph>
void WriteComponents(const Graph& g,
		ComponentSplitter<typename Graph::VertexId> &inner_splitter,
		const AbstractFilter<vector<typename Graph::VertexId>> &checker,
		const string& file_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler, const string& graph_name =
				"my_graph") {
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	WriteComponents<Graph>(g, splitter, file_name, colorer, labeler,
			graph_name);
}

template<class Graph>
void WriteComponent(const GraphComponent<Graph>& gc,
		const string& file_name, const GraphColorer<Graph> &colorer,
		const GraphLabeler<Graph>& labeler, const string& graph_name =
				"my_graph") {
	PrecountedComponentSplitter<typename Graph::VertexId> splitter(gc.v_begin(), gc.v_end());
	WriteComponents(gc.g(), splitter, file_name, colorer, labeler, graph_name);
}

template<class Graph>
void WriteComponents(const Graph& g, size_t split_edge_length,
		const string& file_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler, const string& graph_name =
				"my_graph") {
	ReliableSplitter<Graph> splitter(g, 60, split_edge_length);
	ComponentSizeFilter<Graph> filter(g, split_edge_length, 2, 300);
	WriteComponents<Graph>(g, splitter, filter, file_name, colorer, labeler,
			graph_name);
}

//todo alert!!! magic constants!!!
//todo refactoring of params needed
template<class Graph>
void WriteComponentsAlongPath(const Graph& g, const AbstractFilter<vector<typename Graph::VertexId>> &checker,
		const GraphLabeler<Graph>& labeler, const string& file_name,
		size_t split_edge_length, size_t component_vertex_number,
		const MappingPath<typename Graph::EdgeId>& path,
		const GraphColorer<Graph>& colorer) {
	//	LongEdgesSplitter<Graph> inner_splitter(g, split_edge_length);
	//	ReliableSplitterAlongGenome(g, 60, split_edge_length, MappingPath<EdgeId> genome_path)
	ReliableSplitterAlongPath<Graph> splitter(g, component_vertex_number/*100*//*60*/, split_edge_length, path);
	WriteComponents<Graph>(g, splitter, checker, file_name,
			colorer, labeler);
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g,
		const GraphLabeler<Graph>& labeler, const string& file_name,
		size_t split_edge_length, size_t component_vertex_number,
		const MappingPath<typename Graph::EdgeId>& path,
		const GraphColorer<Graph>& colorer) {
	//	LongEdgesSplitter<Graph> inner_splitter(g, split_edge_length);
	//	ReliableSplitterAlongGenome(g, 60, split_edge_length, MappingPath<EdgeId> genome_path)
	ReliableSplitterAlongPath<Graph> splitter(g, component_vertex_number/*100*//*60*/, split_edge_length, path);
	WriteComponents<Graph>(g, splitter, file_name,
			colorer, labeler);
}

template<class Graph>
void WriteComponentsAlongPath(const Graph& g,
		const GraphLabeler<Graph>& labeler, const string& file_name,
		size_t split_edge_length, size_t component_vertex_number,
		const MappingPath<typename Graph::EdgeId>& path,
		Path<typename Graph::EdgeId> color1 = Path<typename Graph::EdgeId>(),
		Path<typename Graph::EdgeId> color2 = Path<typename Graph::EdgeId>(),
		bool colour_path = false) {
//	Path<typename Graph::EdgeId> simple_path1 = color1.simple_path();
//	Path<typename Graph::EdgeId> simple_path2 = color2.simple_path();
	PathColorer<Graph> path_colorer(g, /*simple_path1*/color1, /*simple_path2*/
	color2);
	auto coloring = path_colorer.ColorPath();
	if (colour_path) {
		for (size_t i = 0; i < path.size(); i++) {
			coloring.insert(make_pair(path[i].first, "green"));
		}
	}
	//	LongEdgesSplitter<Graph> inner_splitter(g, split_edge_length);
	//	ReliableSplitterAlongGenome(g, 60, split_edge_length, MappingPath<EdgeId> genome_path)
	auto_ptr<GraphColorer<Graph>> colorer = DefaultColorer(g, coloring);
	WriteComponentsAlongPath(g, labeler, file_name, split_edge_length, component_vertex_number, path, *colorer);
}

//todo alert!!! magic constants!!!
template<class Graph>
void WriteComponentsAlongGenome(
		const Graph& g,
		const GraphLabeler<Graph>& labeler,
		const string& file_name,
		size_t split_edge_length,
		MappingPath<typename Graph::EdgeId> color1 = MappingPath<
				typename Graph::EdgeId>(),
		MappingPath<typename Graph::EdgeId> color2 = MappingPath<
				typename Graph::EdgeId>()) {
	WriteComponentsAlongPath<Graph>(g, labeler, file_name, split_edge_length, 60,
			color1, color1.simple_path(), color2.simple_path());
}

template<class Graph>
void WriteComponentsAroundEdge(const Graph& g, typename Graph::EdgeId e,
		const string& file_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler, size_t length_bound = 40000, size_t max_size = 30) {
	EdgeNeighborhoodFinder<Graph> splitter(g, e, max_size, length_bound);
	WriteComponents(g, splitter/*, "locality_of_edge_" + ToString(g_.int_id(edge))*/
	, file_name, colorer, labeler);
}

template<class Graph>
void WriteComponentsAroundEdge(const Graph& g, const AbstractFilter<vector<typename Graph::VertexId>> &filter, typename Graph::EdgeId e,
		const string& file_name, const GraphColorer<Graph>& colorer,
		const GraphLabeler<Graph>& labeler) {
	EdgeNeighborhoodFinder<Graph> splitter(g, e, 30, 40000);
	FilteringSplitterWrapper<Graph> filtering_splitter(splitter, filter);
	WriteComponents(g, filtering_splitter/*, "locality_of_edge_" + ToString(g_.int_id(edge))*/
	, file_name, colorer, labeler);
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
