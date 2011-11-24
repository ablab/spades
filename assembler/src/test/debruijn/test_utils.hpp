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

void PrintGraphComponentContainingEdge(const string& file_name, const Graph& g,
		size_t split_edge_length, const IdTrackHandler<Graph>& int_ids) {
	LongEdgesInclusiveSplitter<Graph> inner_splitter(g, split_edge_length);
	ComponentSizeFilter<Graph> checker(g, split_edge_length, 2);
	FilteringSplitterWrapper<Graph> splitter(inner_splitter, checker);
	while (!splitter.Finished()) {
		auto component_vertices = splitter.NextComponent();

		ConjugateDataPrinter<Graph> printer(g, component_vertices.begin(),
				component_vertices.end(), int_ids);
		PrintBasicGraph<Graph>(file_name, printer);
	}
}

}
#endif /* TEST_UTILS_HPP_ */
