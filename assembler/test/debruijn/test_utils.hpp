///*
// * test_utils.hpp
// *
// *  Created on: Apr 10, 2011
// *      Author: sergey
// */
//
//#ifndef TEST_UTILS_HPP_
//#define TEST_UTILS_HPP_
//
//#include "read_generator.hpp"
//#include "edge_graph_constructor.hpp"
//#include "coverage_handler.hpp"
//#include "paired_info.hpp"
//
//namespace de_bruijn_test {
//
//using edge_graph::EdgeGraph;
//
//template <size_t k, class Stream>
//void ConstructGraph(EdgeGraph& graph, Stream& stream) {
//	typedef EdgeGraph::EdgeId EdgeId;
//	typedef de_bruijn::DeBruijnPlus<k + 1, EdgeId> DeBruijn;
//	typedef de_bruijn::PairedInfoIndex<EdgeGraph> PairedIndex;
//
//	SimpleReaderWrapper<Stream> unitedStream(stream);
//	DeBruijn debruijn(unitedStream);
//	EdgeGraph g(k);
//	de_bruijn::EdgeIndex<k+1, EdgeGraph> index(g, debruijn);
//
//
//	stream.reset();
//	edge_graph::EdgeGraphConstructor<k> g_c(debruijn);
//	g_c.ConstructGraph(g, index);
//
//	de_bruijn::CoverageHandler<EdgeGraph> coverageHandler(g);
//	coverageHandler.FillCoverage(unitedStream, index);
//
//	stream.reset();
//	PairedIndex paired_info_index(g, I);
//
//	paired_info_index.FillIndex<K, ReadStream> (index, stream);
//	ClipTips(g, index, genome, "tips_clipped.dot");
//
//	RemoveBulges(g, index, genome, "bulges_removed.dot");
//}
//
//void SimpleTool(istream& is) {
//	size_t size = 100;
//	size_t coverage = 10;
//	size_t gap = 10;
//	ReadGenerator generator(2, size, is, coverage, gap);
//	EdgeGraph g(5);
////	ConstructGraph<5, ReadGenerator>(g, generator);
//}
//
//}
//#endif /* TEST_UTILS_HPP_ */
