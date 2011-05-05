/*
 * coverageCounter.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: sergey
 */

#ifndef COVERAGECOUNTER_HPP_
#define COVERAGECOUNTER_HPP_

using namespace std;

#include <algorithm>
#include "edge_graph.hpp"
#include "utils.hpp"

namespace edge_graph {

template<size_t k, class Graph>
class CoverageCounter {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename de_bruijn::EdgeIndex<k + 1, Graph> Index;
private:
	Graph& g_;
	const de_bruijn::SimpleSequenceMapper<k, Graph> threader_;
	const Index& index_;

	void processRead(Read read) {
		de_bruijn::Path<EdgeId> path = threader_.MapSequence(
				Sequence(read.getSequenceString()));
		if (path.sequence().size() == 0)
			return;
		const vector<EdgeId> &sequence = path.sequence();
		for (typename vector<EdgeId>::const_iterator it = sequence.begin(); it
				!= path.sequence().end(); ++it) {
			g_.IncCoverage(*it, g_.length(*it));
		}
		g_.IncCoverage(sequence[0], -path.start_pos());
		Edge *last = sequence[sequence.size() - 1];
		g_.IncCoverage(last, path.end_pos() - g_.length(last));
	}
public:
	CoverageCounter(Graph& g,
			const Index& index) :
		g_(g), threader_(g, index), index_(index) {
	}

	template<typename FwdIt_>
	void CountCoverage(FwdIt_ begin, FwdIt_ end) {
		for (FwdIt_ r_it = begin; r_it != end; ++r_it) {
			processRead(*r_it);
		}
	}

	template<typename stream>
	void CountCoverage(stream &reader) {
		while (!reader.eof()) {
			Read read;
			reader >> read;
			processRead(read);
		}
	}
};
}

#endif /* COVERAGECOUNTER_HPP_ */
