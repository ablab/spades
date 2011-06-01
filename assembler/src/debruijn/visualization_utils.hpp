#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "utils.hpp"
#include "graphVisualizer.hpp"

namespace debruijn_graph {

using gvis::PairedGraphPrinter;

//template<class Graph>
//class VisHandler: public TraversalHandler<Graph> {
//	typedef typename Graph::VertexId VertexId;
//	typedef typename Graph::EdgeId EdgeId;
//	const Graph& g_;
//	gvis::GraphPrinter<typename Graph::VertexId>& pr_;
//public:
//
//	VisHandler(const Graph& g, gvis::GraphPrinter<VertexId>& pr) :
//		g_(g), pr_(pr) {
//	}
//
//	virtual void HandleVertex(VertexId v) {
//		stringstream ss;
//		ss << v << " comp " << g_.conjugate(v) << "   ";
//		pr_.addVertex(v, ss.str());
//	}
//
//	virtual void HandleEdge(EdgeId e) {
//		stringstream ss;
//		ss << e << " comp " << g_.conjugate(e) << " len "
//				<< g_.EdgeNucls(e).size() << "   ";
//		pr_.addEdge(g_.EdgeStart(e), g_.EdgeEnd(e), ss.str());
//	}
//
//};
//
//template<class Graph>
//class ConjugateVisHandler: public TraversalHandler<Graph> {
//	typedef typename Graph::VertexId VertexId;
//	typedef typename Graph::EdgeId EdgeId;
//	const Graph& g_;
//	gvis::PairedGraphPrinter<VertexId>& pr_;
//	const map<EdgeId, string> color_;
//	string ConstructLabel(EdgeId e) {
//		stringstream ss;
//		//		if (g_.length(e) > 10)
//		ss << g_.length(e);
//		//		else
//		//			ss << g_.length(e) << ":" << e->nucls();
//		ss << "(";
//		ss << ((int) (g_.coverage(e) * 100)) * 0.01;
//		ss << ") id=" << e;
//		return ss.str();
//	}
//public:
//
//	ConjugateVisHandler(const Graph& g, gvis::PairedGraphPrinter<VertexId>& pr,
//			map<EdgeId, string> color) :
//		g_(g), pr_(pr), color_(color) {
//	}
//
//	ConjugateVisHandler(const Graph& g, gvis::PairedGraphPrinter<VertexId>& pr) :
//		g_(g), pr_(pr), color_() {
//	}
//
//	virtual void HandleVertex(VertexId v) {
//		pr_.addVertex(v, "", g_.conjugate(v), "");
//	}
//
//	virtual void HandleEdge(EdgeId e) {
//		VertexId v1 = g_.EdgeStart(e);
//		VertexId v2 = g_.EdgeEnd(e);
//		auto col = color_.find(e);
//		if (col == color_.end())
//			pr_.addEdge(make_pair(v1, g_.conjugate(v1)),
//					make_pair(v2, g_.conjugate(v2)), ConstructLabel(e));
//		else
//			pr_.addEdge(make_pair(v1, g_.conjugate(v1)),
//					make_pair(v2, g_.conjugate(v2)), ConstructLabel(e),
//					col->second);
//	}
//
//};

template<class Graph>
class GraphVisualizer {
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
	gvis::GraphPrinter<VertexId>& gp_;
public:
	GraphVisualizer(const Graph& g, gvis::GraphPrinter<VertexId>& gp) :
		g_(g), gp_(gp) {

	}

	virtual ~GraphVisualizer() {

	}

	virtual void Visualize() = 0;

};

template<class Graph>
class PartialGraphVisualizer {
	typedef typename Graph::VertexId VertexId;
	const Graph& g_;
	gvis::GraphPrinter<VertexId>& gp_;
public:
	PartialGraphVisualizer(const Graph& g, gvis::GraphPrinter<VertexId>& gp) :
		g_(g), gp_(gp) {
	}

	virtual ~PartialGraphVisualizer() {

	}

	virtual void Visualize(const vector<VertexId>& vertices) = 0;

};

//template<class Graph>
//class SimpleGraphVisualizer: public GraphVisualizer {
//	gvis::GraphPrinter<typename Graph::VertexId>& gp_;
//	Graph graph_;
//public:
//	SimpleGraphVisualizer(gvis::GraphPrinter<typename Graph::VertexId>& gp,
//			Graph &graph) :
//		gp_(gp), graph_(graph) {
//	}
//
//	virtual void Visualize() {
//		VisHandler < Graph > h(graph_, gp_);
//		DFS<Graph> (graph_).Traverse(&h);
//		//		gp_.output();
//	}
//};
//
//template<class Graph>
//class ConjugateGraphVisualizer: public GraphVisualizer {
//	gvis::PairedGraphPrinter<typename Graph::VertexId>& gp_;
//	Graph &graph_;
//public:
//	ConjugateGraphVisualizer(
//			gvis::PairedGraphPrinter<typename Graph::VertexId>& gp, Graph graph) :
//		gp_(gp), graph_(graph) {
//	}
//
//	virtual void Visualize() {
//		ConjugateVisHandler < Graph > h(graph_, gp_);
//		DFS<Graph> (graph_).Traverse(&h);
//		gp_.output();
//	}
//};

template<class Graph>
class AdapterGraphVisualizer: GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef GraphVisualizer<Graph> super;
	PartialGraphVisualizer<Graph>& partial_visualizer_;
public:
	AdapterGraphVisualizer(PartialGraphVisualizer<Graph>& partial_visualizer) :
		partial_visualizer_(partial_visualizer) {

	}

	virtual void Visualize() {
		partial_visualizer_.Visualize(
				vector<VertexId> (super::g_.begin(), super::g_.end()));
	}
};

template<class Graph>
class ColoredGraphVisualizer: GraphVisualizer<Graph> {
	typedef GraphVisualizer<Graph> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const map<EdgeId, size_t> &edge_color_;
	//	Graph &graph_;

public:
	ColoredGraphVisualizer(gvis::PairedGraphPrinter<VertexId>& gp,
			map<EdgeId, size_t> &edge_color, Graph &graph) :
		super(graph, gp), edge_color_(edge_color) {
	}

	virtual void Visualize(const vector<VertexId> &vertices) {
		//TODO
		//		gp_.output();
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

	map<EdgeId, string> ColorPath(const debruijn_graph::Path<EdgeId> path) {
		map<EdgeId, string> colors;
		constructColorMap(colors, graph_, path);
		return colors;
	}
};

//template<class Graph>
//class ColoredPathGraphVisualizer: ColoredGraphVisualizer<Graph> {
//	typedef typename Graph::VertexId VertexId;
//	typedef typename Graph::EdgeId EdgeId;
//	gvis::PairedGraphPrinter<VertexId>& gp_;
//	const map<EdgeId, size_t> &edge_color_;
//	const Path<EdgeId> path_;
//
//public:
//	ColoredPathGraphVisualizer(gvis::PairedGraphPrinter<VertexId>& gp,
//			const Path<EdgeId> path, const Graph& g) : ColoredGraphVisualizer<Graph>(),
//		gp_(gp), path_(path) {
//
//	}
//
//	virtual void Visualize(const vector<VertexId> vertices) {
//		PathColorer colorer(g);
//		map<EdgeId, size_t> colors = colorer.ColorPath(path);
//		visualizer.Visualize(g, g.begin(), g.end());
//	}
//};

template<class Graph>
void WriteToFile(const string& file_name, const string& graph_name,
		const Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId> ()) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<typename Graph::VertexId> gp(
			"simulated_data_graph", filestr);
//	ColoredPathGraphVisualizer<Graph> gv(gp, path);
//	gv.Visualize(g);
//	filestr.close();
//	string simple_file_name(file_name);
//	simple_file_name.insert(simple_file_name.size() - 4, "_simple");
//	filestr.open((simple_file_name).c_str(), fstream::out);
//	gvis::GraphPrinter<typename Graph::VertexId> gpr("simulated_data_graph",
//			filestr);
//	SimpleGraphVisualizer<Graph> sgv(gpr);
//	sgv.Visualize(g);
	filestr.close();
}

template<class Graph>
class SmartVisualizer {
private:
	Graph& graph_;
	Sequence genome_;
	const string name_;
public:
	SmartVisualizer(Graph &graph, Sequence genome, string name) :
		graph_(graph), genome_(genome), name_(name) {
	}

	void Visualize(ostream &out = cout) {
		PairedGraphPrinter<typename Graph::VertexId> printer(name_, &out);

	}
};

}

#endif /* VISUALIZATIONUTILS_HPP_ */
