#ifndef VISUALIZATIONUTILS_HPP_
#define VISUALIZATIONUTILS_HPP_

#include "edge_graph.hpp"
#include "coverage_counter.hpp"

namespace edge_graph {

class ColoredPathGraphVisualizer: public GraphVisualizer {
	gvis::PairedGraphPrinter<VertexId>& gp_;
	const de_bruijn::Path<EdgeId> path_;

	void SetColor(map<EdgeId, string> &color, Edge *edge, string col) {
		map<EdgeId, string>::iterator it = color.find(edge);
		if(it != color.end() && it->second != col) {
			color[edge] = "purple";
		} else
			color[edge] = col;
	}

	void constructColorMap(map<EdgeId, string> &color, const EdgeGraph &g, const de_bruijn::Path<EdgeId> path) {
		for (de_bruijn::Path<EdgeId>::iterator it = path.sequence().begin(); it
				!= path.sequence().end(); ++it) {
			SetColor(color, *it, "red");
			Edge* e = *it;
			EdgeId edge = g.Complement(e);
			SetColor(color, edge, "blue");
		}
	}

public:
	ColoredPathGraphVisualizer(gvis::PairedGraphPrinter<VertexId>& gp,
			const de_bruijn::Path<EdgeId> path) :
		gp_(gp), path_(path) {
	}

	virtual void Visualize(const EdgeGraph& g) {
		map<EdgeId, string> color;
		constructColorMap(color, g, path_);
		ComplementVisHandler h(g, gp_, color);
		de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
		gp_.output();
	}
};

}

#endif /* VISUALIZATIONUTILS_HPP_ */
