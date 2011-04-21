#define K 27
#define DE_BRUIJN_DATA_FOLDER "./data/debruijn/"

/////////////////
//for read generator
//#define SUBSTR_LENGTH 10000
//#define COVERAGE 30
//#define R 35
/////////////////


#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "debruijn_graph_test.hpp"

#include "edge_graph_test.hpp"
#include "edge_graph_tool.hpp"
#include "visualization_utils.hpp"
#include "ifaststream.hpp"

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
//	s += de_bruijn::DeBruijnGraphSuite();
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

void RunEdgeGraphTool() {
	pair<pair<string, string>, int> input = QUAKE_CROPPED_10_5;
//	ireadstream stream1(input.first.first);
//	ireadstream stream2(input.first.second);
	string reads[2] = {input.first.first, input.first.second};
	StrobeReader<2, Read, ireadstream> reader((string *)reads);
	ifaststream genome_stream(ECOLI_FILE);
	string genome;
	genome_stream >> genome >> genome;
	edge_graph::EdgeGraphTool(reader, genome.substr(0, input.second));
	reader.close();
	genome_stream.close();
}

int main() {
//	RunTestSuites();
	RunEdgeGraphTool();
	return 0;
}
