#define SUBSTR_LENGTH 10000
#define COVERAGE 30
#define R 35
#define K 27

#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
//#include "condensed_graph_test.hpp"
//#include "condensed_graph_tool.hpp"
#include "debruijn_graph_test.hpp"
#include "edge_graph_test.hpp"
#include "edge_graph_tool.hpp"
#include "visualization_utils.hpp"
#include "ifaststream.hpp"

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
//	s += DeBruijnGraphSuite();
//	s += condensed_graph::CondensedGraphSuite();
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

void RunEdgeGraphTool() {
	ireadstream stream(QUAKE_CROPPED_10_4_A);
	ifaststream genome_stream(ECOLI_FILE);
	string genome;
	genome_stream >> genome >> genome;
	edge_graph::EdgeGraphTool(stream, genome);
	stream.close();
	genome_stream.close();
}

int main() {
//	RunTestSuites();
	RunEdgeGraphTool();
	return 0;
}
