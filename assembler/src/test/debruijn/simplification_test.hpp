#pragma once

#include <boost/test/unit_test.hpp>
#include "test_utils.hpp"
#include "graph_simplification.hpp"

namespace debruijn_graph {

BOOST_AUTO_TEST_SUITE(graph_simplification_tests)

	static debruijn_config::simplification::bulge_remover standard_config_generation() {
		debruijn_config::simplification::bulge_remover br_config;
		br_config.max_length_div_K = 3;
		br_config.max_coverage = 1000.;
		br_config.max_relative_coverage = 1.2;
		br_config.max_delta = 3;
		br_config.max_relative_delta = 0.1;
		return br_config;
	}

	debruijn_config::simplification::bulge_remover standard_config() {
		static debruijn_config::simplification::bulge_remover br_config = standard_config_generation();
		return br_config;
	}

BOOST_AUTO_TEST_CASE( SimpleBulgeRemovalTest ) {
	Graph g(55);
	IdTrackHandler<Graph> int_ids(g);
	ScanBasicGraph("./src/test/debruijn/graph_fragments/simpliest_bulge/simpliest_bulge", g, int_ids);
	RemoveBulges(g, standard_config());
	BOOST_CHECK_EQUAL(g.size(), 4);
}

BOOST_AUTO_TEST_SUITE_END()
}
