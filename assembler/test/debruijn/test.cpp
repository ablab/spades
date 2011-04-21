#define K 27
#define R 100
#define I 220
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
#include "ireadstream.hpp"

#include <tr1/tuple>

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
//	s += de_bruijn::DeBruijnGraphSuite();
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

void RunEdgeGraphTool() {
	const tr1::tuple<string, string, int> input = QUAKE_CROPPED_10_5;
	const string reads[2] = {tr1::get<0>(input), tr1::get<1>(input)};
	StrobeReader<2, Read, ireadstream> reader(reads);
	ireadstream genome_stream(ECOLI_FILE);
	Read genome;
	genome_stream >> genome;
	edge_graph::EdgeGraphTool(reader,  genome.getSequenceString().substr(0, tr1::get<2>(input)));
	reader.close();
	genome_stream.close();
}

int main() {
//	RunTestSuites();
	RunEdgeGraphTool();
	return 0;
}
