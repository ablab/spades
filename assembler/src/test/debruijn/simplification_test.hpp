#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "graph_simplification.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(graph_simplification_tests)

BOOST_AUTO_TEST_CASE( SimpleBulgeRemovalTest ) {
	Graph g(55);
	IdTrackHandler<Graph> int_ids(g);
	ScanBasicGraph("/home/snurk/git/algorithmic-biology/assembler/src/test/debruijn/graph_fragments/simpliest_bulge", g, int_ids);
	RemoveBulges(g);
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		cout << g.length(*it) << endl;
	}
	BOOST_CHECK_EQUAL(g.size(), 4);
}

BOOST_AUTO_TEST_SUITE_END()
}
