#include "visualization_utils.hpp"

namespace edge_graph {

void WriteToFile(const string& file_name, const string& graph_name,
		const EdgeGraph& g, de_bruijn::Path<EdgeId> path) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<VertexId> gp("simulated_data_graph", filestr);
	ColoredPathGraphVisualizer gv(gp, path);
	gv.Visualize(g);
	filestr.close();
}

void SimpleGraphVisualizer::Visualize(const EdgeGraph& g) {
	VisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}

void ComplementGraphVisualizer::Visualize(const EdgeGraph& g) {
	ComplementVisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}

}
