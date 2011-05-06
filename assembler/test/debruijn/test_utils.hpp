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

namespace de_bruijn_test {

using edge_graph::EdgeGraph;
using de_bruijn::EdgeIndex;
using de_bruijn::PairedInfoIndex;

template <size_t k, class ReadStream>
void ConstructGraph(EdgeGraph& g, PairedInfoIndex<EdgeGraph>& paired_index, ReadStream& read_stream) {
	typedef EdgeGraph::EdgeId EdgeId;
	typedef de_bruijn::DeBruijnPlus<k + 1, EdgeId> DeBruijn;
	typedef SimpleReaderWrapper<ReadStream> UnitedStream;
	UnitedStream unitedStream(read_stream);
	DeBruijn debruijn(unitedStream);
	EdgeIndex<k + 1, EdgeGraph> index(g, debruijn);
	edge_graph::CondenseGraph<k> (debruijn, g, index);
	de_bruijn::CoverageHandler<EdgeGraph> coverage_handler(g);
	edge_graph::FillCoverage<k, UnitedStream> (coverage_handler, unitedStream, index);
	edge_graph::FillPairedIndex<k, ReadStream> (paired_index, read_stream, index);
}

template <size_t k>
void ConstructGraphFromGenome(EdgeGraph& g, PairedInfoIndex<EdgeGraph>& paired_index, const string& genome, size_t read_size) {
	typedef read_generator::ReadGenerator<read_generator::SmoothPositionChooser> Stream;
	size_t coverage = read_size;
	size_t gap = 0;
	Stream raw_stream(2, read_size, genome, coverage, gap);
	typedef RCReaderWrapper<Stream> RCStream;
	RCStream read_stream(raw_stream);
	de_bruijn_test::ConstructGraph<k, RCStream>(g, paired_index, read_stream);
}

}
#endif /* TEST_UTILS_HPP_ */
