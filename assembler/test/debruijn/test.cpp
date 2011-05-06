#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "edge_graph_test.hpp"
#include "test_utils.hpp"

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
	s += edge_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

using edge_graph::EdgeGraph;

int main() {
//	RunTestSuites();
	string genome = "AAAAAAAAAAAAAAAAAAAAA";
	EdgeGraph g(5);
	de_bruijn::PairedInfoIndex<EdgeGraph> paired_index(g);
	de_bruijn_test::ConstructGraphFromGenome<5>(g, paired_index, genome, 10);
	return 0;
}
