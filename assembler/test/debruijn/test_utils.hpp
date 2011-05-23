/*
 * test_utils.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include "read_generator.hpp"
#include "launch.hpp"

namespace debruijn_graph {

//using edge_graph::EdgeGraph;
template <size_t k>
void ConstructGraphFromGenome(EdgeGraph& g, EdgeIndex<k + 1, EdgeGraph>& index, CoverageHandler<EdgeGraph>& coverage_handler, PairedInfoIndex<EdgeGraph>& paired_index, const string& genome, size_t read_size) {
	typedef read_generator::ReadGenerator<read_generator::SmoothPositionChooser> Stream;
	size_t coverage = 2*read_size;
	size_t gap = 0;
	Stream raw_stream(2, read_size, genome, coverage, gap);
	typedef RCReaderWrapper<Stream, PairedRead> RCStream;
	RCStream read_stream(raw_stream);
	ConstructGraphWithPairedInfo<k, RCStream>(g, index, coverage_handler, paired_index, read_stream);
}

}
#endif /* TEST_UTILS_HPP_ */
