/*
 * test_utils.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

namespace de_bruijn_test {

template <size_t k, class Stream>
void ConstructGraph(EdgeGraph& graph, Stream& stream) {
	SimpleReaderWrapper<ReadStream> unitedStream(stream);
	DeBruijn debruijn(unitedStream);
	EdgeGraph g(k);
	de_bruijn::DeBruijnPlus<k + 1, EdgeId> &index = debruijn;
	de_bruijn::EdgeHashRenewer<k + 1, EdgeGraph> index_handler(g, index);

	de_bruijn::CoverageHandler<EdgeGraph> coverageHandler(g);

	stream.reset();
	CondenseGraph<ReadStream> (debruijn, g, index, stream, genome);

	stream.reset();
	PairedIndex paired_info_index(g, I);

	FillPairedIndex(paired_info_index, stream, index);
	ClipTips(g, index, genome, "tips_clipped.dot");

	RemoveBulges(g, index, genome, "bulges_removed.dot");
}

}
#endif /* TEST_UTILS_HPP_ */
