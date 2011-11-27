#define BOOST_TEST_MODULE debruijn_test

#include "graphio.hpp"
#include <iostream>
#include "logging.hpp"
#include "test_utils.hpp"

//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"

DECL_PROJECT_LOGGER("dt")

namespace debruijn_graph {

//BOOST_AUTO_TEST_CASE( GenerateGraphFragment ) {
//	std::string input_path = "/home/snurk/git/algorithmic-biology/assembler/data/debruijn/";
//	std::string output_path = "/home/snurk/git/algorithmic-biology/assembler/src/test/debruijn_test/graph_fragments";
//	size_t split_threshold = 100;
//	int int_edge_id = 100;
//	conj_graph_pack gp((Sequence()));
//	PairedInfoIndex<Graph> clustered_index(gp.g);
//	ScanWithClusteredIndex(input_path, gp, clustered_index);
//	PrintGraphComponentContainingEdge(output_path, gp.g,
//			split_threshold, gp.int_ids, int_edge_id);
//}

}

