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

#include "condensed_graph.hpp"
#include <algorithm>

namespace condensed_graph {

class Path {
	vector<Vertex*> vertices_;
	size_t start_pos_;
	size_t end_pos_;

public:

	Path(vector<Vertex*> vertices, size_t start_pos, size_t end_pos) :
		vertices_(vertices), start_pos_(start_pos), end_pos_(end_pos) {

	}

	size_t start_pos() const {
		return start_pos_;
	}

	size_t end_pos() const {
		return end_pos_;
	}

	const vector<Vertex*>& vertices() const {
		return vertices_;
	}
};

template<size_t k>
class SimpleReadThreader {
	const CondensedGraph& g_;
	const SimpleIndex<k>& index_;
public:
	SimpleReadThreader(const CondensedGraph& g, const SimpleIndex<k>& index) :
		g_(g), index_(index) {
	}

	template <typename T>
	Path ThreadRead(const T& read) {
		vector<Vertex*> passed;
		Seq<k> kmer = read.start();
		assert(index_.contains(kmer));
		pair<Vertex*, size_t> graph_start_pos = index_.get(kmer);
		Vertex* last_passed = graph_start_pos.first;
		passed.push_back(last_passed);
		for (size_t i = k; i < read.size() - 1; ++i) {
			kmer = k << read[i];
			assert(index_.contains(kmer));
			Vertex* v = index_.get(k).first;
			if (v != last_passed) {
				passed.push_back(v);
				last_passed = v;
			}
		}
		kmer = read.end();
		pair<Vertex*, size_t> graph_end_pos = index_.get(kmer);
		if (graph_end_pos.first != last_passed) {
			passed.push_back(graph_end_pos.first);
//			last_passed = graph_end_pos.first;
		}
		return Path(passed, graph_start_pos.second, graph_end_pos.second + k);
	}
};

template<size_t k>
class CoverageCounter {
	const CondensedGraph& g_;
	const SimpleIndex<k>& index_;
public:
	CoverageCounter(const CondensedGraph& g, const SimpleIndex<k>& index) :
		g_(g), index_(index) {
	}

	template <typename T, typename FwdIt_>
	void CountCoverage(FwdIt_ begin, FwdIt_ end) {
		SimpleReadThreader<k> threader(g_, index_);

		for (FwdIt_ r_it = begin; r_it != end; ++r_it) {
			Path path = threader.ThreadRead(*r_it);
			vector<Vertex*>::const_iterator v_it = path.vertices().start();
			Vertex* prev = *v_it;
			++v_it;
			while (v_it != path.vertices().end()) {
				Vertex* curr = *v_it;
				char kth_nucl = curr->nucls()[k];
				prev -> inc_coverage(kth_nucl);
				prev = curr;
				++v_it;
			}
		}
	}
};

}
#endif /* COVERAGECOUNTER_HPP_ */
