#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "graph_simplification.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(graph_simplification_tests)

BOOST_AUTO_TEST_CASE( SimpleBulgeRemovalTest ) {
	conj_graph_pack gp;
	ScanGraphPack("", gp);
	RemoveBulges(gp.g);
}

BOOST_AUTO_TEST_SUITE_END()
}
