#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
//#include "condensedGraphTest.hpp"
//#include "condensed_graph_tool.hpp"
#include "debruijnGraphTest.hpp"
#include "edgeGraphTest.hpp"
#include "tip_clipper.hpp"

void runSuite() {
	cute::suite s;
	//TODO add your test here
	s += DeBruijnGraphSuite();
	//	 s += CondensedGraphSuite();
	s += EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suite");
}

void checkClipTippingCompilation() {
	EdgeGraph graph(11);
	TipComparator comparator(graph);
	TipClipper<TipComparator> clipper(graph, comparator, 3, 2.);
	clipper.ClipTips();
}

int main() {
	runSuite();
	//	 SimulatedMistakesTool();
	checkClipTippingCompilation();
	return 0;
}
