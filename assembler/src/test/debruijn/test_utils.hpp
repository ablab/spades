/*
 * test_utils.hpp
 *
 *  Created on: Apr 10, 2011
 *      Author: sergey
 */

#ifndef TEST_UTILS_HPP_
#define TEST_UTILS_HPP_

#include "read/read_generator.hpp"
#include "launch.hpp"

namespace debruijn_graph {

template <size_t k>
void ConstructGraphFromGenome(Graph& g, EdgeIndex<k + 1, Graph>& index/*, CoverageHandler<DeBruijnGraph>& coverage_handler*/
		, PairedInfoIndex<Graph>& paired_index, const string& genome, size_t read_size) {
	typedef read_generator::ReadGenerator<read_generator::SmoothPositionChooser> Stream;
	size_t coverage = 2*read_size;
	size_t gap = 0;
	Stream raw_stream(2, read_size, genome, coverage, gap);
	typedef io::RCReaderWrapper<io::PairedRead> RCStream;
	RCStream read_stream(raw_stream);
	ConstructGraphWithPairedInfo<k, RCStream>(g, index/*, coverage_handler*/, paired_index, read_stream);
}

}
#endif /* TEST_UTILS_HPP_ */
