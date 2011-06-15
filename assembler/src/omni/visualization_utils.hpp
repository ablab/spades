#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "graph_printer.hpp"
#include "graph_labeler.hpp"
#include "omni_utils.hpp"
#include "stack"
#include "queue"

namespace gvis {

//DECL_LOGGER("omg.gvis")

using gvis::PairedGraphPrinter;
using omnigraph::SmartEdgeIterator;
using omnigraph::Path;

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
			omnigraph::GraphLabeler<Graph>& gl) :
		super(g), gp_(gp), gl_(gl) {
	}

	virtual void Visualize() {
		gp_.open();
		for (auto it = super::g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			gp_.AddVertex(*it, gl_.label(*it));
		}
		for (auto it = super::g_.SmartEdgeBegin(); !it.IsEnd(); ++it) {
			gp_.AddEdge(super::g_.EdgeStart(*it), super::g_.EdgeEnd(*it), gl_.label(*it));
		}
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
			for (auto e_it = edges.begin(); e_it != edges.end(); ++e_it) {
				VertexId edge_end = super::g_.EdgeEnd(*e_it);
				if (vertex_set.count(edge_end) > 0) {
					super::gp_.AddEdge(*v_it, edge_end, gl_.label(*e_it), EdgeColor(*e_it));
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
	Path<EdgeId> path_;

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) {
		auto it = color.find(edge);
		if (it != color.end() && it->second != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void ConstructColorMap(map<EdgeId, string> &color) {
		for (auto it = path_.sequence().begin(); it != path_.sequence().end(); ++it) {
			SetColor(color, *it, "red");
			EdgeId e = *it;
			EdgeId edge = graph_.conjugate(e);
			SetColor(color, edge, "blue");
		}
	}

	void ConstructBlackEdgesSet(set<EdgeId> &result) {
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator) {
			result.insert(*iterator);
		}
		for (auto iterator = path_.sequence().begin(); iterator
				!= path_.sequence().end(); ++iterator) {
			result.erase(*iterator);
			result.erase(graph_.conjugate(*iterator));
		}
	}

public:
	PathColorer(Graph &graph, Path<EdgeId> path) :
		graph_(graph), path_(path) {
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

template<class Element>
class GraphSplitter {
public:
	virtual vector<Element> NextComponent() = 0;

	virtual bool Finished() = 0;

	virtual ~GraphSplitter() {

	}
};

template<class Graph>
class ErrorComponentSplitter: public GraphSplitter<typename Graph::VertexId> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph &graph_;
	set<EdgeId> black_edges_;
	SmartEdgeIterator<Graph> iterator_;
	set<VertexId> visited_;

public:
	ErrorComponentSplitter(Graph &graph, const set<EdgeId> &black_edges) :
		graph_(graph), black_edges_(black_edges),
				iterator_(graph.SmartEdgeBegin()) {
	}

	virtual ~ErrorComponentSplitter() {
	}

	size_t FindComponent(VertexId start_vertex, set<VertexId> &component) {
		std::stack<VertexId> st;
		st.push(start_vertex);
		size_t result = 0;
		while (!st.empty()) {
			VertexId next = st.top();
			st.pop();
			if (visited_.find(next) == visited_.end()) {
				component.insert(next);
				visited_.insert(next);
				vector < EdgeId > outgoing = graph_.OutgoingEdges(next);
				for (size_t i = 0; i < outgoing.size(); i++) {
					if (black_edges_.find(outgoing[i]) != black_edges_.end()) {
						//						black_edges_.erase(outgoing[i]);
						result += graph_.length(outgoing[i]);
						st.push(graph_.EdgeEnd(outgoing[i]));
					}
				}
				vector < EdgeId > incoming = graph_.IncomingEdges(next);
				//				cout << outgoing.size() << " " << incoming.size() << " "
				//						<< next << endl;
				for (size_t i = 0; i < incoming.size(); i++) {
					if (black_edges_.find(incoming[i]) != black_edges_.end()) {
						//						black_edges_.erase(incoming[i]);
						result += graph_.length(incoming[i]);
						st.push(graph_.EdgeStart(incoming[i]));
					}
				}
			}
		}
		return result;
	}

	void ProcessOutgoing(std::priority_queue<pair<size_t, VertexId> > &q,
			set<VertexId> &was, VertexId next, size_t next_length, size_t bound) {
		vector < EdgeId > outgoing = graph_.OutgoingEdges(next);
		for (size_t i = 0; i < outgoing.size(); i++) {
			VertexId v = graph_.EdgeEnd(outgoing[i]);
			if (was.find(v) == was.end()) {
				size_t size = next_length + graph_.length(outgoing[i]);
				q.push(make_pair(size, v));
			}
		}
	}

	void ProcessIncoming(std::priority_queue<pair<size_t, VertexId> > &q,
			set<VertexId> &was, VertexId next, size_t next_length, size_t bound) {
		vector < EdgeId > incoming = graph_.IncomingEdges(next);
		for (size_t i = 0; i < incoming.size(); i++) {
			VertexId v = graph_.EdgeStart(incoming[i]);
			if (was.find(v) == was.end()) {
				size_t size = next_length + graph_.length(incoming[i]);
				q.push(make_pair(size, graph_.EdgeStart(incoming[i])));
			}
		}
	}

	void Dijkstra(set<VertexId> &component, size_t bound) {
		set < VertexId > was;
		std::priority_queue < pair<size_t, VertexId> > q;
		for (auto iterator = component.begin(); iterator != component.end(); ++iterator) {
			q.push(make_pair(0, *iterator));
		}
		while (!q.empty()) {
			auto next_pair = q.top();
			q.pop();
			VertexId next = next_pair.second;
			size_t next_length = next_pair.first;
			if (was.find(next) != was.end()) {
				continue;
			}
			component.insert(next_pair.second);
			was.insert(next_pair.second);
			if (next_length <= bound) {
				ProcessOutgoing(q, was, next, next_length, bound);
				ProcessIncoming(q, was, next, next_length, bound);
			}
		}
	}

	virtual vector<VertexId> NextComponent() {
		if (Finished()) {
			assert(false);
			return vector<VertexId> ();
		}
		EdgeId next = *iterator_;
		++iterator_;
		set < VertexId > component;
		//		cout << "oppa" << endl;
		size_t component_size = FindComponent(graph_.EdgeEnd(next), component);
		//		cout << component.size() << endl;
		Dijkstra(component, component_size);
		//		cout << component.size() << endl;
		return vector<VertexId> (component.begin(), component.end());
	}

	virtual bool Finished() {
		while (!iterator_.IsEnd()) {
			if (black_edges_.find(*iterator_) != black_edges_.end()
					&& visited_.find(graph_.EdgeEnd(*iterator_))
							== visited_.end()) {
				return false;
			}
			++iterator_;
		}
		return true;
	}

};

template<class Graph>
void WriteSimple(const string& file_name, const string& graph_name, Graph& g, GraphLabeler<Graph> labeler) {
	fstream filestr;
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	gvis::DotGraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(g, gpr, labeler);
	sgv.Visualize();
	filestr.close();
}

template<class Graph>
void WritePaired(const string& file_name, const string& graph_name, Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId> ()) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
	PathColorer<Graph> path_colorer(g, path);
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	omnigraph::StrGraphLabeler<Graph> gl(g);
	ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
	AdapterGraphVisualizer<Graph> result_vis(g, gv);
	result_vis.Visualize();
	filestr.close();
}

template<class Graph>
string ConstructComponentName(string file_name, size_t cnt) {
	stringstream ss;
	ss << "_error_" << cnt;
	string res = file_name;
	res.insert(res.length() - 4, ss.str());
//	cout << res << endl;
	return res;
}

template<class Graph>
void WriteErrors(const string& file_name, const string& graph_name, Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId> ()) {
	PathColorer<Graph> path_colorer(g, path);
	set<typename Graph::EdgeId> black = path_colorer.BlackEdges();
	omnigraph::StrGraphLabeler<Graph> gl(g);
	ErrorComponentSplitter<Graph> splitter(g, black);
	size_t cnt = 0;
	map<typename Graph::EdgeId, string> coloring = path_colorer.ColorPath();
	//	cout << "oppa" << endl;
	while (!splitter.Finished() && cnt < 100) {
		fstream filestr;
		filestr.open(ConstructComponentName<Graph> (file_name, cnt).c_str(),
				fstream::out);
		gvis::DotPairedGraphPrinter<Graph> gp(g, graph_name, filestr);
		ColoredGraphVisualizer<Graph> gv(g, gp, gl, coloring);
		auto component = splitter.NextComponent();
		gp.open();
		gv.Visualize(component);
		gp.close();
		//		cout << cnt << endl;
		cnt++;
	}
}

template<class Graph>
void WriteToFile(const string& file_name, const string& graph_name, Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId> ()) {
	WritePaired(file_name, graph_name, g, path);
	WriteSimple(file_name, graph_name, g);
	WriteErrors(file_name, graph_name, g, path);
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
