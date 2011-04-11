#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "condensed_graph_test.hpp"
#include "condensed_graph_tool.hpp"
#include "debruijn_graph_test.hpp"
#include "edge_graph_test.hpp"
#include "edge_graph_tool.hpp"
#include "tip_clipper.hpp"

void runSuite() {
	cute::suite s;
	//TODO add your test here
	s += DeBruijnGraphSuite();
	s += condensed_graph::CondensedGraphSuite();
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suite");
}

void checkClipTippingCompilation() {
	using namespace de_bruijn;
	using namespace edge_graph;
	EdgeGraph graph(11);
	TipComparator comparator(graph);
	TipClipper<TipComparator> clipper(comparator, 3, 2.);
	clipper.ClipTips(graph);
}

int main() {
//	runSuite();
	edge_graph::SimulatedMistakesTool();
	//	 SimulatedMistakesTool();
//	checkClipTippingCompilation();
	return 0;
}
