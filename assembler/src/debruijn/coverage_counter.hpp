/*
 * coverageCounter.hpp
 *
 *  Created on: Mar 18, 2011
 *      Author: sergey
 */

#ifndef COVERAGECOUNTER_HPP_
#define COVERAGECOUNTER_HPP_

using namespace std;

#include "condensed_graph.hpp"
#include <algorithm>
#include "edge_graph.hpp"
#include "utils.hpp"

namespace edge_graph {

template<size_t k, class Graph>
class SimpleReadThreader {
public:
	typedef typename Graph::EdgeId EdgeId;
private:
	const Graph& g_;
	const de_bruijn::SimpleIndex<k + 1, EdgeId>& index_;

	void processKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition) const {
		if (index_.contains(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			endPosition = position.second;
			if (passed.empty()) {
				startPosition = position.second;
			}
			if (passed.empty() || passed[passed.size() - 1] != position.first)
				passed.push_back(position.first);
		}
	}
public:
	SimpleReadThreader(const Graph& g, const de_bruijn::SimpleIndex<k + 1, EdgeId>& index) :
		g_(g), index_(index) {
	}

	de_bruijn::Path<EdgeId> ThreadRead(const Sequence& read) const {
		vector<EdgeId> passed;
		if(read.size() <= k)
			return de_bruijn::Path<EdgeId>();
		Seq<k + 1> kmer = read.start<k + 1> ();
		size_t startPosition = -1;
		size_t endPosition = -1;
		processKmer(kmer, passed, startPosition, endPosition);
		for (size_t i = k; i < read.size(); ++i) {
			processKmer(kmer, passed, startPosition, endPosition);
			kmer = kmer << read[i];
		}
		return de_bruijn::Path<EdgeId> (passed, startPosition, endPosition);
	}
};

template<size_t k, class Graph>
class CoverageCounter {
public:
	typedef typename Graph::EdgeId EdgeId;
private:
	Graph& g_;
	const SimpleReadThreader<k, Graph> threader_;
	const de_bruijn::SimpleIndex<k + 1, EdgeId>& index_;

	void processRead(Read read) {
		de_bruijn::Path<EdgeId> path = threader_.ThreadRead(
				Sequence(read.getSequenceString()));
		if (path.sequence().size() == 0)
			return;
		const vector<EdgeId> &sequence = path.sequence();
		for (typename vector<EdgeId>::const_iterator it = sequence.begin(); it
				!= path.sequence().end(); ++it) {
			g_.inc_coverage(*it, g_.length(*it));
		}
		g_.inc_coverage(sequence[0], -path.start_pos());
		Edge *last = sequence[sequence.size() - 1];
		g_.inc_coverage(last, path.end_pos() - g_.length(last));
	}
public:
	CoverageCounter(Graph& g, const de_bruijn::SimpleIndex<k + 1, EdgeId>& index) :
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

//namespace condensed_graph {
//
//typedef edge_graph::Path Path;
////class Path {
////	vector<Vertex*> vertices_;
////	size_t start_pos_;
////	size_t end_pos_;
////
////public:
////
////	Path(vector<Vertex*> vertices, size_t start_pos, size_t end_pos) :
////		vertices_(vertices), start_pos_(start_pos), end_pos_(end_pos) {
////
////	}
////
////	size_t start_pos() const {
////		return start_pos_;
////	}
////
////	size_t end_pos() const {
////		return end_pos_;
////	}
////
////	const vector<Vertex*>& vertices() const {
////		return vertices_;
////	}
////};
//
//template<size_t k>
//class SimpleReadThreader {
//	const CondensedGraph& g_;
//	const SimpleIndex<k>& index_;
//public:
//	SimpleReadThreader(const CondensedGraph& g, const SimpleIndex<k>& index) :
//		g_(g), index_(index) {
//	}
//
//	template<typename T>
//	Path ThreadRead(const T& read) {
//		vector<Vertex*> passed;
//		Seq<k + 1> kmer = read.start();
//		assert(index_.contains(kmer));
//		pair<Vertex*, size_t> graph_start_pos = index_.get(kmer);
//		Vertex* last_passed = graph_start_pos.first;
//		passed.push_back(last_passed);
//		for (size_t i = k + 1; i < read.size() - 1; ++i) {
//			kmer = kmer << read[i];
//
//			Vertex* v = index_.get(k).first;
//			if (v != last_passed) {
//				passed.push_back(v);
//				last_passed = v;
//			}
//		}
//		kmer = read.end();
//		pair<Vertex*, size_t> graph_end_pos = index_.get(kmer);
//		if (graph_end_pos.first != last_passed) {
//			passed.push_back(graph_end_pos.first);
//			//			last_passed = graph_end_pos.first;
//		}
//		return Path(passed, graph_start_pos.second, graph_end_pos.second + k);
//	}
//};
//
//template<size_t k>
//class CoverageCounter {
//	const CondensedGraph& g_;
//	const SimpleIndex<k>& index_;
//public:
//	CoverageCounter(const CondensedGraph& g, const SimpleIndex<k>& index) :
//		g_(g), index_(index) {
//	}
//
//	template<typename T, typename FwdIt_>
//	void CountCoverage(FwdIt_ begin, FwdIt_ end) {
//		SimpleReadThreader<k> threader(g_, index_);
//
//		for (FwdIt_ r_it = begin; r_it != end; ++r_it) {
//			Path path = threader.ThreadRead(*r_it);
//			vector<Vertex*>::const_iterator v_it = path.vertices().start();
//			Vertex* prev = *v_it;
//			++v_it;
//			while (v_it != path.vertices().end()) {
//				Vertex* curr = *v_it;
//				char kth_nucl = curr->nucls()[k];
//				prev -> inc_coverage(kth_nucl);
//				prev = curr;
//				++v_it;
//			}
//		}
//	}
//};
//
//}
#endif /* COVERAGECOUNTER_HPP_ */
