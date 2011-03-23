/*
 * coverageCounter.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: sergey
 */

#ifndef COVERAGECOUNTER_HPP_
#define COVERAGECOUNTER_HPP_

using namespace std;
using namespace condensed_graph;

class CoverageCounter {
	CondensedGraph *g_;
	SimpleHashTable *h_t_;
public:
	CoverageCounter(CondensedGraph *g, SimpleHashTable *h_t) : g_(g), h_t_(h_t) {

	}

	void CountCoverage() {

	}
};

#endif /* COVERAGECOUNTER_HPP_ */
