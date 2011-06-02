#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "utils.hpp"
#include "graphVisualizer.hpp"

namespace debruijn_graph {

template <class Graph>
class VisHandler: public TraversalHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	gvis::GraphPrinter<typename Graph::VertexId>& pr_;
public:

	VisHandler(const Graph& g, gvis::GraphPrinter<VertexId>& pr) :
		g_(g), pr_(pr) {
	}

	virtual void HandleVertex(VertexId v) {
		stringstream ss;
		ss << v<< " comp "<< g_.conjugate(v)<<"   ";
		pr_.addVertex(v, ss.str());
	}

	virtual void HandleEdge(EdgeId e) {
		stringstream ss;
		ss << e<<" comp "<<g_.conjugate(e)<<" len "<<g_.EdgeNucls(e).size()<<"   ";
		pr_.addEdge(g_.EdgeStart(e), g_.EdgeEnd(e), ss.str());
	}

};

template <class Graph>
class ConjugateVisHandler: public TraversalHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;
	gvis::PairedGraphPrinter<VertexId>& pr_;
	const map<EdgeId, string> color_;
	string ConstructLabel(EdgeId e) {
		stringstream ss;
//		if (g_.length(e) > 10)
			ss << g_.length(e);
//		else
//			ss << g_.length(e) << ":" << e->nucls();
		ss << "(";
		ss << ((int) (g_.coverage(e) * 100)) * 0.01;
		ss << ") id="<<e;
		return ss.str();
	}
public:

	ConjugateVisHandler(const Graph& g,
			gvis::PairedGraphPrinter<VertexId>& pr, map<EdgeId, string> color) :
		g_(g), pr_(pr), color_(color) {
	}

	ConjugateVisHandler(const Graph& g,
			gvis::PairedGraphPrinter<VertexId>& pr) :
		g_(g), pr_(pr), color_() {
	}

	virtual void HandleVertex(VertexId v) {
		pr_.addVertex(v, "", g_.conjugate(v), "");
	}

	virtual void HandleEdge(EdgeId e) {
		VertexId v1 = g_.EdgeStart(e);
		VertexId v2 = g_.EdgeEnd(e);
		auto col = color_.find(e);
		if (col == color_.end())
			pr_.addEdge(make_pair(v1, g_.conjugate(v1)),
					make_pair(v2, g_.conjugate(v2)), ConstructLabel(e));
		else
			pr_.addEdge(make_pair(v1, g_.conjugate(v1)),
					make_pair(v2, g_.conjugate(v2)), ConstructLabel(e),
					col->second);
	}

};

template <class Graph>
class GraphVisualizer {
public:
	virtual void Visualize(const Graph& g) = 0;
};

template <class Graph>
class SimpleGraphVisualizer: public GraphVisualizer<Graph> {
	gvis::GraphPrinter<typename Graph::VertexId>& gp_;
public:
	SimpleGraphVisualizer(gvis::GraphPrinter<typename Graph::VertexId>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const Graph& g) {
		VisHandler<Graph> h(g, gp_);
		DFS<Graph>(g).Traverse(&h);
		gp_.output();
	}
};

template <class Graph>
class ConjugateGraphVisualizer: public GraphVisualizer<Graph> {
	gvis::PairedGraphPrinter<typename Graph::VertexId>& gp_;
public:
	ConjugateGraphVisualizer(gvis::PairedGraphPrinter<typename Graph::VertexId>& gp) :
		gp_(gp) {
	}

	virtual void Visualize(const Graph& g) {
		ConjugateVisHandler<Graph> h(g, gp_);
		DFS<Graph>(g).Traverse(&h);
		gp_.output();
	}
};

template <class Graph>
class ColoredPathGraphVisualizer: public GraphVisualizer<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	gvis::PairedGraphPrinter<VertexId>& gp_;
	const debruijn_graph::Path<EdgeId> path_;

	void SetColor(map<EdgeId, string> &color, EdgeId edge, string col) {
		auto it = color.find(edge);
		if(it != color.end() && it->second != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void constructColorMap(map<EdgeId, string> &color, const Graph &g, const Path<EdgeId> path) {
		for (auto it = path.sequence().begin(); it
				!= path.sequence().end(); ++it) {
			SetColor(color, *it, "red");
			EdgeId e = *it;
			EdgeId edge = g.conjugate(e);
			SetColor(color, edge, "blue");
		}
	}

public:
	ColoredPathGraphVisualizer(gvis::PairedGraphPrinter<VertexId>& gp,
			const debruijn_graph::Path<EdgeId> path) :
		gp_(gp), path_(path) {
	}

	virtual void Visualize(const Graph& g) {
		map<EdgeId, string> color;
		constructColorMap(color, g, path_);
		ConjugateVisHandler<Graph> h(g, gp_, color);
		debruijn_graph::DFS<Graph>(g).Traverse(&h);
		gp_.output();
	}
};

template <class Graph>
void WriteToFile(const string& file_name, const string& graph_name,
		const Graph& g,
		Path<typename Graph::EdgeId> path = Path<typename Graph::EdgeId>()) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<typename Graph::VertexId> gp(graph_name, filestr);
	ColoredPathGraphVisualizer<Graph> gv(gp, path);
	gv.Visualize(g);
	filestr.close();
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size()-4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	gvis::GraphPrinter<typename Graph::VertexId> gpr(graph_name, filestr);
	SimpleGraphVisualizer<Graph> sgv(gpr);
	sgv.Visualize(g);
	filestr.close();
}

}

#endif /* VISUALIZATIONUTILS_HPP_ */
