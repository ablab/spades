#include "visualization_utils.hpp"


void debruijn_graph::WriteToFile(const string& file_name, const string& graph_name,
		const EdgeGraph& g, de_bruijn::Path<EdgeId> path) {
	fstream filestr;
	filestr.open(file_name.c_str(), fstream::out);
	gvis::PairedGraphPrinter<VertexId> gp("simulated_data_graph", filestr);
	ColoredPathGraphVisualizer gv(gp, path);
	gv.Visualize(g);
	filestr.close();
	string simple_file_name(file_name);
	simple_file_name.insert(simple_file_name.size()-4, "_simple");
	filestr.open((simple_file_name).c_str(), fstream::out);
	gvis::GraphPrinter<VertexId> gpr("simulated_data_graph", filestr);
	SimpleGraphVisualizer sgv(gpr);
	sgv.Visualize(g);
	filestr.close();


}

void debruijn_graph::SimpleGraphVisualizer::Visualize(const EdgeGraph& g) {
	VisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}

void debruijn_graph::ConjugateGraphVisualizer::Visualize(const EdgeGraph& g) {
	ConjugateVisHandler h(g, gp_);
	de_bruijn::DFS<EdgeGraph>(g).Traverse(&h);
	gp_.output();
}
