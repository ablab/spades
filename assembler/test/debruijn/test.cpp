#include "cute.h"
#include "ide_listener.h"
#include "cute_runner.h"
#include "debruijn_graph_test.hpp"
#include "test_utils.hpp"

void RunTestSuites() {
	cute::suite s;
	//TODO add your test here
	s += debruijn_graph::EdgeGraphSuite();
	cute::ide_listener lis;
	cute::makeRunner(lis)(s, "De Bruijn Project Test Suites");
}

int main() {
	RunTestSuites();
//	string genome = "AAAAAAAAAAAAAAAAAAAAA";
//	EdgeGraph g(5);
//	de_bruijn::EdgeIndex<5 + 1, EdgeGraph> index(g);
//	de_bruijn::CoverageHandler<EdgeGraph> coverage_handler(g);
//	de_bruijn::PairedInfoIndex<EdgeGraph> paired_index(g);
//	edge_graph::ConstructGraphFromGenome<5>(g, index, coverage_handler, paired_index, genome, 10);
	return 0;
}
