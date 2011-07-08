#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "graph_printer.hpp"
#include "graph_labeler.hpp"
#include "omni_utils.hpp"
#include "stack"
#include "queue"
#include "dijkstra.hpp"
#include "observable_graph.hpp"
#include "splitters.hpp"

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
protected:
	Graph& g_;
	//	const GraphLabeler<Graph> &labeler_;
	//	gvis::GraphPrinter<VertexId>& gp_;
public:
	GraphVisualizer(Graph& g/*, gvis::GraphPrinter<VertexId>& gp*/) :
		g_(g) {

	}

	virtual ~GraphVisualizer() {

	}

	virtual void Visualize() = 0;

};

template<class Graph>
class PartialGraphVisualizer {
	typedef typename Graph::VertexId VertexId;
protected:
	Graph& g_;
	gvis::GraphPrinter<VertexId>& gp_;
public:
	PartialGraphVisualizer(Graph& g, gvis::GraphPrinter<VertexId>& gp) :
		g_(g), gp_(gp) {
	}

	void open() {
		gp_.open();
	}

	void close() {
		gp_.close();
	}

	virtual ~PartialGraphVisualizer() {

	}

	virtual void Visualize(const vector<VertexId>& vertices) = 0;

};

template<class Graph>
class SimpleGraphVisualizer: public GraphVisualizer<Graph> {
	typedef GraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	gvis::GraphPrinter<VertexId>& gp_;
	const omnigraph::GraphLabeler<Graph>& gl_;
public:
	SimpleGraphVisualizer(Graph& g, gvis::GraphPrinter<VertexId>& gp,
			const omnigraph::GraphLabeler<Graph>& gl) :
		super(g), gp_(gp), gl_(gl) {
	}

	virtual void Visualize() {
		gp_.open();
		TRACE("Visualize started");
		for (auto it = super::g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			gp_.AddVertex(*it, gl_.label(*it));
		}
		TRACE("Vertices printed");
		for (auto it = super::g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			gp_.AddEdge(super::g_.EdgeStart(*it), super::g_.EdgeEnd(*it),
					gl_.label(*it));
		}
		TRACE("Edges printed");
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
	AdapterGraphVisualizer(Graph& g,
			PartialGraphVisualizer<Graph>& partial_visualizer) :
		super(g), partial_visualizer_(partial_visualizer) {

	}

	virtual void Visualize() {
		partial_visualizer_.open();
		partial_visualizer_.Visualize(
				vector<VertexId> (super::g_.begin(), super::g_.end()));
		partial_visualizer_.close();
	}
};

template<class Graph>
class ColoredGraphVisualizer: public PartialGraphVisualizer<Graph> {
	typedef PartialGraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const omnigraph::GraphLabeler<Graph>& gl_;
	const map<EdgeId, string>& edge_colors_;
	const string default_color_;
	const string border_vertex_color_;

	bool IsBorder(VertexId v, const set<VertexId>& vertices) {
		const vector<EdgeId> outgoing_edges = super::g_.OutgoingEdges(v);
		const vector<EdgeId> incoming_edges = super::g_.IncomingEdges(v);
		set < EdgeId > adjacent_edges;
		adjacent_edges.insert(outgoing_edges.begin(), outgoing_edges.end());
		adjacent_edges.insert(incoming_edges.begin(), incoming_edges.end());
		for (auto e_it = adjacent_edges.begin(); e_it != adjacent_edges.end(); ++e_it) {
			if (vertices.count(super::g_.EdgeStart(*e_it)) == 0
					|| vertices.count(super::g_.EdgeEnd(*e_it)) == 0) {
				return true;
			}
		}
		return false;
	}

	string EdgeColor(EdgeId e) {
		string edge_color = default_color_;
		auto edge_colors_it = edge_colors_.find(e);
		if (edge_colors_it != edge_colors_.end()) {
			edge_color = (*edge_colors_it).second;
		}
		return edge_color;
	}

public:
	ColoredGraphVisualizer(Graph& g, gvis::GraphPrinter<VertexId>& gp,
			const omnigraph::GraphLabeler<Graph>& gl,
			const map<EdgeId, string>& edge_colors,
			const string& default_color = "",
			const string& border_vertex_color = "yellow") :
		super(g, gp), gl_(gl), edge_colors_(edge_colors),
				default_color_(default_color),
				border_vertex_color_(border_vertex_color) {
	}

	virtual void Visualize(const vector<VertexId>& vertices) {
		set < VertexId > vertex_set(vertices.begin(), vertices.end());
		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			string vertex_color =
					IsBorder(*v_it, vertex_set) ? border_vertex_color_
							: "white";
			super::gp_.AddVertex(*v_it, gl_.label(*v_it), vertex_color);
		}
		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			const vector<EdgeId> edges = super::g_.OutgoingEdges(*v_it);
			TRACE("working with vertex " << *v_it);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = super::g_.EdgeEnd(*e_it);
				TRACE(super::g_.coverage(*e_it)<<" " << super::g_.length(*e_it));
				if (vertex_set.count(edge_end) > 0) {
					super::gp_.AddEdge(*v_it, edge_end, gl_.label(*e_it),
							EdgeColor(*e_it));
					TRACE("Edge added");
				}
			}
		}
	}

};

template<class Graph>
class PathColorer {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	Graph &graph_;
	Path<EdgeId> path1_;
	Path<EdgeId> path2_;

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) {
		auto it = color.find(edge);
		if (it != color.end() && it->second != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void ConstructColorMap(map<EdgeId, string> &color) {
		for (auto it = path1_.sequence().begin(); it != path1_.sequence().end(); ++it) {
			SetColor(color, *it, "red");
		}
		for (auto it = path2_.sequence().begin(); it != path2_.sequence().end(); ++it) {
			SetColor(color, *it, "blue");
		}
	}

	void ConstructBlackEdgesSet(set<EdgeId> &result) {
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator) {
			result.insert(*iterator);
		}
		for (auto iterator = path1_.sequence().begin(); iterator
				!= path1_.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
		for (auto iterator = path2_.sequence().begin(); iterator
				!= path2_.sequence().end(); ++iterator) {
			result.erase(*iterator);
		}
	}

public:
	PathColorer(Graph &graph, Path<EdgeId> path1, Path<EdgeId> path2) :
		graph_(graph), path1_(path1), path2_(path2) {
	}

	map<EdgeId, string> ColorPath() {
		map<EdgeId, string> colors;
		ConstructColorMap(colors);
		return colors;
	}

	set<EdgeId> BlackEdges() {
		set<EdgeId> result;
		ConstructBlackEdgesSet(result);
		return result;
	}
};

template<class Graph>
void WriteToDotFile(const string& file_name, const string& graph_name,
		Graph& g,
		const GraphLabeler<Graph>& labeler = EmptyGraphLabeler<Graph> ()) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(g, gpr, labeler);
	sgv.Visualize();
	filestr.close();
}

template<class Graph>
void WriteSimple(const string& file_name, const string& graph_name, Graph& g,
		const GraphLabeler<Graph>& labeler) {
	DEBUG("Writing simple graph "<<file_name);
	fstream filestr;
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	gvis::DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	DEBUG("DotGraphPrinter created");
	SimpleGraphVisualizer<Graph> sgv(g, gpr, labeler);
	DEBUG("SimpleGraphVisualizer created");
	sgv.Visualize();
	DEBUG("sgv.Visualize() ok");
	filestr.close();
	DEBUG("WriteSimple ok");
}

template<class Graph>
void WriteSimple(const string& file_name, const string& graph_name, Graph& g,
		Path<typename Graph::EdgeId> path1, Path<typename Graph::EdgeId> path2) {
	fstream filestr;
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
	filestr.open(simple_file_name.c_str(), fstream::out);
	gvis::DotGraphPrinter<typename Graph::VertexId> gp(graph_name, filestr);
	PathColorer<Graph> path_colorer(g, path1, path2);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	omnigraph::StrGraphLabeler<Graph> gl(g);
	ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

template<class Graph>
void WritePaired(
		const string& file_name,
		const string& graph_name,
		Graph& g,
		Path<typename Graph::EdgeId> path1/* = Path<typename Graph::EdgeId> ()*/,
		Path<typename Graph::EdgeId> path2/* = Path<typename Graph::EdgeId> ()*/) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
	PathColorer<Graph> path_colorer(g, path1, path2);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	omnigraph::StrGraphLabeler<Graph> gl(g);
	ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

template<class Graph>
class AbstractVisualizerFactory {
public:
	virtual PartialGraphVisualizer<Graph> *GetVisualizerInstance(
			gvis::GraphPrinter<typename Graph::VertexId> &gp) = 0;
	virtual gvis::GraphPrinter<typename Graph::VertexId> *GetPrinterInstance(
			const string &graph_name, ostream &os) = 0;
	virtual ~AbstractVisualizerFactory() {
	}
};

template<class Graph>
class ColoredVisualizerFactory: public AbstractVisualizerFactory<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Graph &graph_;
	const GraphLabeler<Graph> &labeler_;
	const map<EdgeId, string> &coloring_;
public:
	ColoredVisualizerFactory(Graph& graph, const GraphLabeler<Graph> &labeler,
			const map<EdgeId, string> &coloring) :
		graph_(graph), labeler_(labeler), coloring_(coloring) {
	}

	virtual gvis::GraphPrinter<typename Graph::VertexId> *GetPrinterInstance(
			const string &graph_name, ostream &os) {
		return new gvis::DotPairedGraphPrinter<Graph>(graph_, graph_name, os);
	}

	virtual PartialGraphVisualizer<Graph> *GetVisualizerInstance(
			gvis::GraphPrinter<typename Graph::VertexId> &gp) {
		return new ColoredGraphVisualizer<Graph> (graph_, gp, labeler_,
				coloring_);
	}

	virtual ~ColoredVisualizerFactory() {
	}
};

template<class Graph>
class ComponentGraphVisualizer: public GraphVisualizer<Graph> {
private:
	AbstractVisualizerFactory<Graph> &factory_;
	GraphSplitter<typename Graph::VertexId> &splitter_;
	const string &file_name_;
	const string &graph_name_;
	size_t max_parts_number_;

	string ConstructComponentName(string file_name, size_t cnt) {
		stringstream ss;
		ss << cnt;
		string res = file_name;
		res.insert(res.length() - 4, ss.str());
		return res;
	}

public:
	ComponentGraphVisualizer(Graph &graph,
			AbstractVisualizerFactory<Graph> &factory,
			GraphSplitter<typename Graph::VertexId> &splitter,
			const string &file_name, const string &graph_name,
			size_t max_parts_number = 100) :
		GraphVisualizer<Graph> (graph), factory_(factory), splitter_(splitter),
				file_name_(file_name), graph_name_(graph_name),
				max_parts_number_(max_parts_number) {
	}

	virtual ~ComponentGraphVisualizer() {
	}

	virtual void Visualize() {
		size_t cnt = 1;
		while (!splitter_.Finished() && cnt <= max_parts_number_) {
			string component_name =
					ConstructComponentName(file_name_, cnt).c_str();
			ofstream os;
			os.open(component_name);
			gvis::GraphPrinter<typename Graph::VertexId> * gp =
					factory_.GetPrinterInstance(graph_name_, os);
			auto visualizer = factory_.GetVisualizerInstance(*gp);
			auto component = splitter_.NextComponent();
			visualizer->open();
			visualizer->Visualize(component);
			visualizer->close();
			os.close();
			delete visualizer;
			delete gp;
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
void WriteErrors(const string& file_name, const string& graph_name, Graph& g,
		Path<typename Graph::EdgeId> path1 = Path<typename Graph::EdgeId> (),
		Path<typename Graph::EdgeId> path2 = Path<typename Graph::EdgeId> ()) {
	PathColorer<Graph> path_colorer(g, path1, path2);
	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
	omnigraph::StrGraphLabeler<Graph> labeler(g);
	ErrorComponentSplitter<Graph> splitter(g, black);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	ColoredVisualizerFactory<Graph> factory(g, labeler, coloring);
	string error_file_name = InsertComponentName<Graph>(file_name, "error");
	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, error_file_name,
			graph_name, 200);
	gv.Visualize();
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

template<class Graph>
void WriteComponents(const string& file_name, const string& graph_name, Graph& g, size_t split_edge_length,
		Path<typename Graph::EdgeId> path1 = Path<typename Graph::EdgeId> (),
		Path<typename Graph::EdgeId> path2 = Path<typename Graph::EdgeId> ()) {
	omnigraph::StrGraphLabeler<Graph> labeler(g);
	PathColorer<Graph> path_colorer(g, path1, path2);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	LongEdgesSplitter<Graph> splitter(g, split_edge_length);
	ColoredVisualizerFactory<Graph> factory(g, labeler, coloring);
	cout << file_name << endl;
	ComponentGraphVisualizer<Graph> gv(g, factory, splitter, file_name,
			graph_name, 200);
	gv.Visualize();
}

template<class Graph>
void WriteToFile(
		const string& file_name,
		const string& graph_name,
		Graph& g,
		Path<typename Graph::EdgeId> path1/* = Path<typename Graph::EdgeId> ()*/,
		Path<typename Graph::EdgeId> path2/* = Path<typename Graph::EdgeId> ()*/) {
	WritePaired(file_name, graph_name, g, path1, path2);
	WriteSimple(file_name, graph_name, g, path1, path2);
	WriteErrors(file_name, graph_name, g, path1, path2);
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
