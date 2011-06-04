#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "graphVisualizer.hpp"
#include "omni_utils.hpp"

namespace omnigraph {

using gvis::PairedGraphPrinter;

template<class Graph>
class GraphVisualizer {
	typedef typename Graph::VertexId VertexId;
protected:
	Graph& g_;
//	gvis::GraphPrinter<VertexId>& gp_;
public:
	GraphVisualizer(Graph& g/*, gvis::GraphPrinter<VertexId>& gp*/) :
		g_(g)/*, gp_(gp)*/ {

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
public:
	SimpleGraphVisualizer(Graph& g, gvis::GraphPrinter<VertexId>& gp)
	: super(g), gp_(gp)
	{}

	virtual void Visualize() {
		gp_.open();
		for (auto it = super::g_.SmartVertexBegin(); !it.isEnd(); ++it) {
			gp_.AddVertex(*it);
		}
		for (auto it = super::g_.SmartEdgeBegin(); !it.isEnd(); ++it) {
			gp_.AddEdge(super::g_.EdgeStart(*it), super::g_.EdgeEnd(*it));
		}
		gp_.close();
	}
};

template<class Graph>
class AdapterGraphVisualizer: public GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef GraphVisualizer<Graph> super;
	PartialGraphVisualizer<Graph>& partial_visualizer_;
public:
	AdapterGraphVisualizer(Graph& g, PartialGraphVisualizer<Graph>& partial_visualizer) :
		super(g), partial_visualizer_(partial_visualizer) {

	}

	virtual void Visualize() {
		partial_visualizer_.open();
		partial_visualizer_.Visualize(vector<VertexId> (super::g_.begin(), super::g_.end()));
		partial_visualizer_.close();
	}
};

template<class Graph>
class ColoredGraphVisualizer: public PartialGraphVisualizer<Graph> {
	typedef PartialGraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const map<EdgeId, string>& edge_colors_;
	const string& default_color_;
	const string& border_vertex_color_;

	bool IsBorder(VertexId v, const set<VertexId>& vertices) {
		const vector<EdgeId> outgoing_edges = super::g_.OutgoingEdges(v);
		const vector<EdgeId> incoming_edges = super::g_.IncomingEdges(v);
		set<EdgeId> adjacent_edges;
		adjacent_edges.insert(outgoing_edges.begin(), outgoing_edges.end());
		adjacent_edges.insert(incoming_edges.begin(), incoming_edges.end());
		for (auto e_it = adjacent_edges.begin(); e_it != adjacent_edges.end(); ++e_it) {
			if (vertices.count(super::g_.EdgeStart(*e_it)) == 0 || vertices.count(super::g_.EdgeEnd(*e_it)) == 0) {
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
			const map<EdgeId, string>& edge_colors, const string& default_color = "", const string& border_vertex_color = "grey") :
		super(g, gp), edge_colors_(edge_colors), default_color_(default_color), border_vertex_color_(border_vertex_color) {
	}

	virtual void Visualize(const vector<VertexId>& vertices) {
		set<VertexId> vertex_set(vertices.begin(), vertices.end());

		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			super::gp_.AddVertex(*v_it, "", IsBorder(*v_it, vertex_set) ? border_vertex_color_:"white");
		}

		for (auto v_it = vertex_set.begin(); v_it != vertex_set.end(); ++v_it) {
			const vector<EdgeId> edges = super::g_.OutgoingEdges(*v_it);
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = super::g_.EdgeEnd(*e_it);
				if (vertex_set.count(edge_end) > 0) {
					super::gp_.AddEdge(*v_it, edge_end, EdgeColor(*e_it));
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

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) {
		auto it = color.find(edge);
		if (it != color.end() && it->second != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void constructColorMap(map<EdgeId, string> &color, const Graph &g,
			const Path<EdgeId> path) {
		for (auto it = path.sequence().begin(); it != path.sequence().end(); ++it) {
			SetColor(color, *it, "red");
			EdgeId e = *it;
			EdgeId edge = g.conjugate(e);
			SetColor(color, edge, "blue");
		}
	}

public:
	PathColorer(Graph &graph) :
		graph_(graph) {
	}

	map<EdgeId, string> ColorPath(const Path<EdgeId> path) {
		map<EdgeId, string> colors;
		constructColorMap(colors, graph_, path);
		return colors;
	}
};

template<class Graph>
void WriteToFile(const string& file_name, const string& graph_name,
		/*const */Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId> ()) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
	PathColorer<Graph> path_colorer(g);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath(path);
	ColoredGraphVisualizer<Graph> gv(g, gp, coloring);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size()-4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	gvis::DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(g, gpr);
	sgv.Visualize();
	filestr.close();
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
