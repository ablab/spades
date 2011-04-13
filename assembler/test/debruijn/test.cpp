#define SUBSTR_LENGTH 10000
#define COVERAGE 30
#define R 35
#define K 15

#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "condensed_graph_test.hpp"
#include "condensed_graph_tool.hpp"
#include "debruijn_graph_test.hpp"
#include "edge_graph_test.hpp"
#include "edge_graph_tool.hpp"
#include "tip_clipper.hpp"
#include "coverage_counter.hpp"

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
//	s += DeBruijnGraphSuite();
//	s += condensed_graph::CondensedGraphSuite();
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

//void checkClipTippingCompilation() {
//	using namespace de_bruijn;
//	using namespace edge_graph;
//	EdgeGraph graph(11);
//	TipComparator comparator(graph);
//	TipClipper<TipComparator> clipper(comparator, 3, 2.);
//	clipper.ClipTips(graph);
//}
<<<<<<< HEAD:assembler/test/debruijn/test.cpp

int main() {
	runSuite();
//	ireadstream stream(QUAKE_CROPPED_10_4_A);
//	edge_graph::ConstructionTool(stream);
//	stream.close();
=======

void RunEdgeGraphConstructionTool() {
	ireadstream stream(ECOLI_FILE);
	edge_graph::ConstructionTool(stream);
	stream.close();
}

int main() {
//	RunTestSuites();
	RunEdgeGraphConstructionTool();
>>>>>>> 7fd1100e920f1d86fca7dcd79637eaa633e8973c:assembler/test/debruijn/test.cpp
	return 0;
}
