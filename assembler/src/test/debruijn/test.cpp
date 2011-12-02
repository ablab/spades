#define BOOST_TEST_MODULE debruijn_test

#include "graphio.hpp"
#include <iostream>
#include "logging.hpp"
#include "test_utils.hpp"


//headers with tests
#include "debruijn_graph_test.hpp"
#include "simplification_test.hpp"
//#include "pair_info_test.hpp"

DECL_PROJECT_LOGGER("dt")

namespace debruijn_graph {

//BOOST_AUTO_TEST_CASE( GenerateGraphFragment ) {
//	std::string input_path = "/home/snurk/git/algorithmic-biology/assembler/data/debruijn/QUAKE_CROPPED_400K/K55/latest/tip_clipping_0/graph";
//	std::string output_path = "/home/snurk/git/algorithmic-biology/assembler/src/test/debruijn/graph_fragments/simpliest_bulge";
//	size_t split_threshold = 1000;
//	int int_edge_id = 1867;
//	graph_pack<ConjugateDeBruijnGraph, 55> gp;
//	ScanGraphPack(input_path, gp);
//	//prints only basic graph structure
//	PrintGraphComponentContainingEdge(output_path, gp.g,
//			split_threshold, gp.int_ids, int_edge_id);
//
//	//long way to write to dot file
//	Graph g(55);
//	IdTrackHandler<Graph> int_ids(g);
//	ScanBasicGraph(output_path, g, int_ids);
//	total_labeler_graph_struct graph_struct(g, &int_ids, (const EdgesPositionHandler<Graph>*)0);
//	total_labeler tot_lab(&graph_struct);
//	WriteToDotFile(g,
//			tot_lab, output_path + ".dot",
//			"mygraph", Path<EdgeId>(), Path<EdgeId>());
//}

}

