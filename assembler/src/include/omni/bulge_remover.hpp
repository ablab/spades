//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * bulge_remover.hpp
 *
 *  Created on: Apr 13, 2011
 *      Author: sergey
 */

#ifndef BULGE_REMOVER_HPP_
#define BULGE_REMOVER_HPP_

#include <cmath>
#include "omni_utils.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"
#include "path_processor.hpp"

namespace omnigraph {

template<class Graph>
struct SimplePathCondition {
	typedef typename Graph::EdgeId EdgeId;
	const Graph& g_;

	SimplePathCondition(const Graph& g) :
			g_(g) {

	}

	bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
		if (edge == g_.conjugate(edge))
			return false;
		for (size_t i = 0; i < path.size(); ++i)
			if (edge == path[i] || edge == g_.conjugate(path[i]))
				return false;
		for (size_t i = 0; i < path.size(); ++i) {
			if (path[i] == g_.conjugate(path[i])) {
				return false;
			}
			for (size_t j = i + 1; j < path.size(); ++j)
				if (path[i] == path[j] || path[i] == g_.conjugate(path[j]))
					return false;
		}
		return true;
	}
};

//template<class Graph>
//struct TrivialCondition {
//	typedef typename Graph::EdgeId EdgeId;
//
//	bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
//		for (size_t i = 0; i < path.size(); ++i)
//			for (size_t j = i + 1; j < path.size(); ++j)
//				if (path[i] == path[j])
//					return false;
//		return true;
//	}
//};

template<class Graph>
bool TrivialCondition(typename Graph::EdgeId edge,
		const vector<typename Graph::EdgeId>& path) {
	typedef typename Graph::EdgeId EdgeId;
	for (size_t i = 0; i < path.size(); ++i)
		for (size_t j = i + 1; j < path.size(); ++j)
			if (path[i] == path[j])
				return false;
	return true;
}

template<class Graph>
class MostCoveredAlternativePathChooser: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	EdgeId forbidden_edge_;
	double max_coverage_;
	vector<EdgeId> most_covered_path_;

	double PathAvgCoverage(const vector<EdgeId>& path) {
		double unnormalized_coverage = 0;
		size_t path_length = 0;
		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId edge = path[i];
			size_t length = g_.length(edge);
			path_length += length;
			unnormalized_coverage += g_.coverage(edge) * length;
		}
		return unnormalized_coverage / path_length;
	}

public:

	MostCoveredAlternativePathChooser(Graph& g, EdgeId edge) :
			g_(g), forbidden_edge_(edge), max_coverage_(-1.0) {

	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		double path_cov = PathAvgCoverage(path);
		for (size_t i = 0; i < path.size(); i++) {
			if (path[i] == forbidden_edge_)
				return;
		}
		if (path_cov > max_coverage_) {
			max_coverage_ = path_cov;
			most_covered_path_ = path;
		}
	}

	double max_coverage() {
		return max_coverage_;
	}

	const vector<EdgeId>& most_covered_path() {
		return most_covered_path_;
	}
};

/**
 * This class removes simple bulges from given graph with the following algorithm: it iterates through all edges of
 * the graph and for each edge checks if this edge is likely to be a simple bulge
 * if edge is judged to be one it is removed.
 */
template<class Graph>
class BulgeRemover {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	typedef boost::function<bool(EdgeId edge, const vector<EdgeId>& path)> BulgeCallbackF;

	/**
	 * Create BulgeRemover with specified parameters.
	 */
	BulgeRemover(Graph& g, size_t max_length, double max_coverage,
			double max_relative_coverage, double max_delta,
			double max_relative_delta, BulgeCallbackF bulge_condition,
			BulgeCallbackF opt_callback = 0,
			boost::function<void(EdgeId)> removal_handler = 0) :
			g_(g), max_length_(max_length), max_coverage_(max_coverage), max_relative_coverage_(
					max_relative_coverage), max_delta_(max_delta), max_relative_delta_(
					max_relative_delta), bulge_condition_(bulge_condition), opt_callback_(
					opt_callback), removal_handler_(removal_handler) {
	}

	~BulgeRemover() {
//		cout << "Time in counter path_processor_time: "
//				<< human_readable_time(path_processor_time.time()) << " number of calls " << path_processor_time.counts() << endl;
//
//		cout << "Time in counter callback_time: "
//				<< human_readable_time(callback_time.time()) << " number of calls " << callback_time.counts() << endl;
//
//		cout << "Time in counter process_bulge_time: "
//				<< human_readable_time(process_bulge_time.time()) << " number of calls " << process_bulge_time.counts() << endl;
//
//		cout << "Time in counter split_time: "
//				<< human_readable_time(split_time.time()) << " number of calls " << split_time.counts() << endl;
//
//		cout << "Time in counter glue_time: "
//				<< human_readable_time(glue_time.time()) << " number of calls " << glue_time.counts() << endl;
//
//		cout << "Time in counter compress_time: "
//				<< human_readable_time(compress_time.time()) << " number of calls " << compress_time.counts() << endl;

	}

	void RemoveBulges();

private:
	Graph& g_;
	size_t max_length_;
	double max_coverage_;
	double max_relative_coverage_;
	double max_delta_;
	double max_relative_delta_;
	BulgeCallbackF bulge_condition_;
	BulgeCallbackF opt_callback_;
	boost::function<void(EdgeId)> removal_handler_;

	//timers
//	avg_perf_counter path_processor_time;
//	avg_perf_counter callback_time;
//
//	avg_perf_counter process_bulge_time;
//	avg_perf_counter split_time;
//	avg_perf_counter glue_time;
//
//	avg_perf_counter compress_time;

	bool PossibleBulgeEdge(EdgeId e);

	size_t PathLength(const vector<EdgeId>& path);

	/**
	 * Returns most covered path from start to the end such that its length doesn't exceed length_left.
	 * Returns pair of empty vector and -1 if no such path could be found.
	 * Edges are returned in reverse order!
	 */
	//	pair<vector<EdgeId> , int> BestPath(VertexId start,
	//			VertexId end, int length_left);
	/**
	 * Checks if alternative path is simple (doesn't contain conjugate edges, edge e or conjugate(e))
	 * and its average coverage is greater than max_relative_coverage_ * g.coverage(e)
	 */
	bool BulgeCondition(EdgeId e, const vector<EdgeId>& path,
			double path_coverage) {
		return math::ge(path_coverage * max_relative_coverage_, g_.coverage(e))
				&& bulge_condition_(e, path);
		//		return path_and_coverage.second > max_relative_coverage * g.kplus_one_mer_coverage(edge);
	}

	void ProcessBulge(EdgeId edge, const vector<EdgeId>& path) {
//		process_bulge_time.start();
//		UniformPositionAligner aligner(PathLength(path) + 1, g_.length(edge) + 1);
		EnsureEndsPositionAligner aligner(PathLength(path), g_.length(edge));
		double prefix_length = 0.;
		vector<size_t> bulge_prefix_lengths;
		for (auto it = path.begin(); it != path.end(); ++it) {
			prefix_length += g_.length(*it);
			bulge_prefix_lengths.push_back(aligner.GetPosition(prefix_length));
		}
		EdgeId edge_to_split = edge;
		size_t prev_length = 0;
		TRACE("Process bulge " << path.size() << " edges");
		for (size_t i = 0; i < path.size(); ++i) {
			if (bulge_prefix_lengths[i] > prev_length) {
				if (bulge_prefix_lengths[i] - prev_length
						!= g_.length(edge_to_split)) {
//					split_time.start();
					pair<EdgeId, EdgeId> split_result = g_.SplitEdge(
							edge_to_split,
							bulge_prefix_lengths[i] - prev_length);
//					split_time.stop();
					edge_to_split = split_result.second;
//					glue_time.start();
					g_.GlueEdges(split_result.first, path[i]);
//					glue_time.stop();
				} else {
//					glue_time.start();
					g_.GlueEdges(edge_to_split, path[i]);
//					glue_time.stop();
				}
			}
			prev_length = bulge_prefix_lengths[i];
		}
//		process_bulge_time.stop();
	}

private:
	DECL_LOGGER("BulgeRemover")
};

template<class Graph>
bool BulgeRemover<Graph>::PossibleBulgeEdge(EdgeId e) {
	return g_.length(e) <= max_length_ && g_.coverage(e) < max_coverage_;
}

template<class Graph>
size_t BulgeRemover<Graph>::PathLength(const vector<EdgeId>& path) {
	size_t length = 0;
	for (size_t i = 0; i < path.size(); ++i) {
		length += g_.length(path[i]);
	}
	return length;
}

template<class Graph>
void BulgeRemover<Graph>::RemoveBulges() {
//	g_.PrintHandlers();
	TRACE("Bulge remove process started");

	size_t it_count = 0;

	DEBUG("RemoveBulges function started");

	CoverageComparator<Graph> comparator(g_);
	for (auto iterator = g_.SmartEdgeBegin(comparator); !iterator.IsEnd();
			++iterator) {

		EdgeId edge = *iterator;

		VERBOSE_POWER_T(++it_count, 1000, "th iteration of bulge processing");

		TRACE(
				"Considering edge " << g_.int_id(edge) << " of length " << g_.length(edge) << " and avg coverage " << g_.coverage(edge));
		TRACE("Is possible bulge " << PossibleBulgeEdge(edge));

		if (PossibleBulgeEdge(edge)) {

			size_t kplus_one_mer_coverage = math::round(
					g_.length(edge) * g_.coverage(edge));

			TRACE(
					"Processing edge " << g_.int_id(edge) << " and coverage " << kplus_one_mer_coverage);

			VertexId start = g_.EdgeStart(edge);
			TRACE("Start " << g_.int_id(start));

			VertexId end = g_.EdgeEnd(edge);
			TRACE("End " << g_.int_id(end));
			size_t delta = std::floor(
					std::max(max_relative_delta_ * g_.length(edge),
							max_delta_));
			MostCoveredAlternativePathChooser<Graph> path_chooser(g_, edge);

//			path_processor_time.start();
			PathProcessor<Graph> path_finder(g_,
					(g_.length(edge) > delta) ? g_.length(edge) - delta : 0,
					g_.length(edge) + delta, start, end, path_chooser);

			path_finder.Process();
//			path_processor_time.stop();

			const vector<EdgeId>& path = path_chooser.most_covered_path();
			double path_coverage = path_chooser.max_coverage();

			TRACE(
					"Best path with coverage " << path_coverage << " is " << PrintPath<Graph>(g_, path));

			//if edge was returned, this condition will fail
			if (BulgeCondition(edge, path, path_coverage)) {

				TRACE("Satisfied condition");

//				callback_time.start();
				if (opt_callback_)
					opt_callback_(edge, path);

				if (removal_handler_)
					removal_handler_(edge);
//				callback_time.stop();

				TRACE("Projecting edge " << g_.int_id(edge));
				ProcessBulge(edge, path);

//				compress_time.start();
				TRACE("Compressing start vertex " << g_.int_id(start))
				g_.CompressVertex(start);

				TRACE("Compressing end vertex " << g_.int_id(end))
				g_.CompressVertex(end);
//				compress_time.stop();

			} else {
				TRACE("Didn't satisfy condition");
			}
		}
		TRACE("-----------------------------------");
	}
}

template<class Graph>
class MostCoveredPathChooser: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	double max_coverage_;
	vector<EdgeId> most_covered_path_;

	double PathAvgCoverage(const vector<EdgeId>& path) {
		double unnormalized_coverage = 0;
		size_t path_length = 0;
		for (size_t i = 0; i < path.size(); ++i) {
			EdgeId edge = path[i];
			size_t length = g_.length(edge);
			path_length += length;
			unnormalized_coverage += g_.coverage(edge) * length;
		}
		return unnormalized_coverage / path_length;
	}

public:

	MostCoveredPathChooser(Graph& g) :
			g_(g), max_coverage_(-1.0) {

	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		double path_cov = PathAvgCoverage(path);
		if (path_cov > max_coverage_) {
			max_coverage_ = path_cov;
			most_covered_path_ = path;
		}
	}

	double max_coverage() {
		return max_coverage_;
	}

	const vector<EdgeId>& most_covered_path() {
		return most_covered_path_;
	}
};

template<class Graph>
class OppositionLicvidator {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	class ComponentFinder {
		Graph& g_;
		size_t max_length_;
		size_t length_difference_;
		VertexId start_v_;
		map<VertexId, pair<size_t, size_t>> processed_;
		set<VertexId> can_be_processed_;
		set<VertexId> neighbourhood_;
		VertexId end_v_;

		bool CanBeProcessed(VertexId v) {
			vector<EdgeId> incoming = g_.IncomingEdges(v);
			TRACE("Check of process possibilities of "<<g_.int_id(v));
			for (auto it = incoming.begin(); it != incoming.end(); ++it) {
				if (processed_.count(g_.EdgeStart(*it)) == 0) {
					TRACE(
							"Blocked by unprocessed or external vertex "<<g_.int_id(g_.EdgeStart(*it))<<" that starts edge "<<g_.int_id(*it));
					return false;
				}
			}
			return true;
		}

		void CountNeighbourhood(VertexId v) {
			vector<EdgeId> outgoing = g_.OutgoingEdges(v);
			for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
				TRACE(
						"Vertex "<<g_.int_id(g_.EdgeEnd(*it)) <<" added to neighbourhood_")
				neighbourhood_.insert(g_.EdgeEnd(*it));
			}
		}

		void CountCanBeProcessedNeighb(VertexId v) {
			vector<EdgeId> outgoing = g_.OutgoingEdges(v);
			for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
				if (CanBeProcessed(g_.EdgeEnd(*it))) {
					TRACE(
							"Vertex "<<g_.int_id(g_.EdgeEnd(*it)) <<" added to can_be_processed_")
					can_be_processed_.insert(g_.EdgeEnd(*it));
				}
			}
		}

		pair<size_t, size_t> DistancesRange(EdgeId e) {
			VertexId start = g_.EdgeStart(e);
			TRACE(
					"Edge "<<g_.int_id(e) <<" of length "<<g_.length(e)<<" with start vertex "<<g_.int_id(start)<<" on distance "<<processed_[g_.EdgeStart(e)]);
			return make_pair(processed_[g_.EdgeStart(e)].first + g_.length(e),
					processed_[g_.EdgeStart(e)].second + g_.length(e));
		}

		void ProcessStartVertex() {
			TRACE("Process vertex "<< g_.int_id(start_v_));
			processed_.insert(make_pair(start_v_, make_pair(0, 0)));
			CountCanBeProcessedNeighb(start_v_);
			CountNeighbourhood(start_v_);
		}

		bool CheckEdgeToStartAbsence(VertexId v) {
			vector<EdgeId> outgoing = g_.OutgoingEdges(v);
			for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
				if (g_.EdgeEnd(*it) == start_v_) {
					return false;
				}
			}
			return true;
		}

		void ProcessVertex(VertexId v) {
			TRACE("Process vertex "<< g_.int_id(v));
			vector<EdgeId> incoming = g_.IncomingEdges(v);
			pair<size_t, size_t> final_range(max_length_, 0);
			for (auto it = incoming.begin(); it != incoming.end(); ++it) {
				pair<size_t, size_t> range = DistancesRange(*it);
				TRACE(
						"Edge "<< g_.int_id(*it) << " provide distance range "<<range);
				if (range.first < final_range.first)
					final_range.first = range.first;
				if (range.second > final_range.second)
					final_range.second = range.second;
			}
			processed_[v] = final_range;
			CountCanBeProcessedNeighb(v);
			CountNeighbourhood(v);
			neighbourhood_.erase(v);
			can_be_processed_.erase(v);
		}

		void ProcessEndVertex(VertexId v) {
			TRACE("Process vertex "<< g_.int_id(v));
			vector<EdgeId> incoming = g_.IncomingEdges(v);
			pair<size_t, size_t> final_range(max_length_, 0);
			for (auto it = incoming.begin(); it != incoming.end(); ++it) {
				pair<size_t, size_t> range = DistancesRange(*it);
				TRACE(
						"Edge "<< g_.int_id(*it) << " provide distance range "<<range);
				if (range.first < final_range.first)
					final_range.first = range.first;
				if (range.second > final_range.second)
					final_range.second = range.second;
			}
			processed_[v] = final_range;
			neighbourhood_.erase(v);
			can_be_processed_.erase(v);
		}

		bool CheckVertexDist(VertexId v) {
			return processed_[v].first < max_length_;
		}

		bool CheckPathLengths() {
			VERIFY(end_v_ != VertexId(NULL));
			return processed_[end_v_].second - processed_[end_v_].first
					< length_difference_;
		}

	public:
		ComponentFinder(Graph& g, size_t max_length, size_t length_difference,
				VertexId start_v) :
				g_(g), max_length_(max_length), length_difference_(
						length_difference), start_v_(start_v), end_v_(
						VertexId(NULL)) {
		}

		bool TryFindComponent() {
			ProcessStartVertex();
			while (neighbourhood_.size() != 1) {
				if (can_be_processed_.empty()) {
					return false;
				} else {
					VertexId v = *(can_be_processed_.begin());
					ProcessVertex(v);
					if (!CheckVertexDist(v) || !CheckEdgeToStartAbsence(v)) {
						return false;
					}
				}
			}
			end_v_ = *(neighbourhood_.begin());
			if (CanBeProcessed(end_v_)) {
				ProcessEndVertex(end_v_);
			} else
				return false;

			return CheckPathLengths();
		}

		const map<VertexId, pair<size_t, size_t>>& processed() const {
			VERIFY(end_v_ != VertexId(NULL));
			return processed_;
		}

		VertexId start_v() const {
			return start_v_;
		}

		VertexId end_v() const {
			return end_v_;
		}

		bool HaveConjugateVertices() const {
			set<VertexId> conjugate_vertices;
			for (auto iter = processed_.begin(); iter != processed_.end();
					++iter) {
				if (conjugate_vertices.find(iter->first)
						== conjugate_vertices.end()) {
					conjugate_vertices.insert(g_.conjugate(iter->first));
				} else {
					return true;
				}
			}
			return false;
		}
	};

	map<VertexId, size_t> AverageDistances(
			const map<VertexId, pair<size_t, size_t>>& ranges) {
		map<VertexId, size_t> answer;
		for (auto it = ranges.begin(); it != ranges.end(); ++it) {
			answer.insert(make_pair(it->first, /*(*/
			it->second.first /*+ it->second.second) / 2*/));
		}
		return answer;
	}

//	EdgeId Project(EdgeId e, EdgeId target, size_t start, size_t end) {
//		EdgeId processed_target = target;
//		bool compress_start = false;
//		bool compress_end = false;
//		if (start > 0) {
//			compress_start = true;
//			pair<EdgeId, EdgeId> split_res = g_.SplitEdge(processed_target,
//					start);
//			processed_target = split_res.second;
//		}
//		if (end < g_.length(target)) {
//			compress_end = true;
//			size_t pos = end > start ? end - start : 1;
//			pair<EdgeId, EdgeId> split_res = g_.SplitEdge(processed_target,
//					pos);
//			processed_target = split_res.first;
//		}
//		EdgeId answer = g_.GlueEdges(e, processed_target);
//		if (compress_start
//				&& g_.CanCompressVertex(g_.EdgeStart(processed_target))) {
//			answer = g_.UnsafeCompressVertex(g_.EdgeStart(processed_target));
//		}
//		if (compress_end
//				&& g_.CanCompressVertex(g_.EdgeEnd(processed_target))) {
//			answer = g_.UnsafeCompressVertex(g_.EdgeEnd(processed_target));
//		}
//		return answer;
//	}

//	EdgeId LicvidateComponent(const map<VertexId, size_t>& component_dist,
//			const vector<EdgeId>& best_path) {
//		for (auto it = component_dist.begin(); it != component_dist.end();
//				++it) {
//			TRACE(
//					"Process vertex "<<g_.int_id(it->first)<<" and distance "<< it->second);
//			VertexId v = it->first;
//			vector<EdgeId> outgoing_edges = g_.OutgoingEdges(v);
//			for (auto e_it = outgoing_edges.begin();
//					e_it != outgoing_edges.end(); ++e_it) {
//				EdgeId e = *e_it;
//				VertexId end_v = g_.EdgeEnd(e);
//				if (e != fake_edge && component_dist.count(end_v) > 0) {
//					TRACE("Project edge "<<g_.int_id(e));
//					fake_edge = Project(e, fake_edge, it->second,
//							component_dist.find(end_v)->second);
//					TRACE(
//							"fake_edge after proj "<<g_.int_id(fake_edge)<< " from "<< g_.int_id(g_.EdgeStart(fake_edge))<<" to "<<g_.int_id(g_.EdgeEnd(fake_edge)));
//					TRACE("Project finished");
//				}
//			}
//		}
//		return fake_edge;
//	}

	set<size_t> AllAvgDist(const map<VertexId, size_t>& avg_dist) {
		set<size_t> answer;
		for (auto it = avg_dist.begin(); it != avg_dist.end(); ++it) {
			answer.insert(it->second);
		}
		return answer;
	}

	bool SinglePath(VertexId start_v, VertexId end_v, size_t min_dist,
			size_t max_dist) {
		PathStorageCallback<Graph> path_storage(g_);
		PathProcessor<Graph> best_path_finder(g_, min_dist, max_dist, start_v,
				end_v, path_storage);
		best_path_finder.Process();
		VERIFY(path_storage.paths().size() > 0);
		return path_storage.paths().size() == 1;
	}

//	EdgeId AddFakeEdge(const vector<EdgeId>& path) {
//		VERIFY(path.size() > 0);
//		vector<const typename Graph::EdgeData*> datas;
//		for (auto it = path.begin(); it != path.end(); ++it) {
//			datas.push_back(&g_.data(*it));
//		}
//		return g_.AddEdge(g_.EdgeStart(path.front()), g_.EdgeEnd(path.back()),
//				g_.master().MergeData(datas));
//	}

	pair<vector<EdgeId>, map<size_t, VertexId>> SplitPath(
			const map<VertexId, size_t>& avg_dist, const vector<EdgeId>& path) {
		VERIFY(!path.empty());
		DEBUG("Splitting path " << g_.str(path));
		vector<EdgeId> split_path;
		map<size_t, VertexId> dist_map;
		set<size_t> all_dist = AllAvgDist(avg_dist);
		for (auto it = path.begin(); it != path.end(); ++it) {
			VertexId start_v = g_.EdgeStart(*it);
			VertexId end_v = g_.EdgeEnd(*it);
			size_t start_dist = avg_dist.find(g_.EdgeStart(*it))->second;
			size_t end_dist = avg_dist.find(g_.EdgeEnd(*it))->second;
			set<size_t> dist_to_split(all_dist.lower_bound(start_dist),
					all_dist.upper_bound(end_dist));
			size_t offset = start_dist;
			EdgeId e = *it;
			for (auto split_it = dist_to_split.begin();
					split_it != dist_to_split.end(); ++split_it) {
				size_t curr = *split_it;
				size_t pos = curr - offset;
				if (pos > 0 && pos < g_.length(e)) {
					DEBUG(
							"Splitting edge " << g_.str(e) << " on position " << pos);
					pair<EdgeId, EdgeId> split_res = g_.SplitEdge(e, pos);
					VertexId inner_v = g_.EdgeEnd(split_res.first);
					DEBUG("Result: edges " << g_.str(split_res.first) << " "
							<< g_.str(split_res.second) << " inner vertex" << inner_v);
					split_path.push_back(split_res.first);
					dist_map[curr] = inner_v;
//					avg_dist[inner_v] = curr;
					e = split_res.second;
					offset = curr;
				}
			}
			split_path.push_back(e);
			dist_map[start_dist] = start_v;
			dist_map[end_dist] = end_v;
		}
		DEBUG("Path splitted");
		return make_pair(split_path, dist_map);
	}

	set<VertexId> KeySet(const map<VertexId, size_t>& avg_dist) {
		set<VertexId> answer;
		for (auto it = avg_dist.begin(); it != avg_dist.end(); ++it) {
			answer.insert(it->first);
		}
		return answer;
	}

	vector<EdgeId> NonPathEdges(const map<VertexId, size_t>& avg_dist,
			const vector<EdgeId>& path) {
		set<VertexId> vertices = KeySet(avg_dist);
		GraphComponent<Graph> component(g_, vertices.begin(), vertices.end());
		set<EdgeId> path_edges(path.begin(), path.end());
		vector<EdgeId> non_path_edges;
		for (auto it = component.e_begin(); it != component.e_end(); ++it) {
			if (path_edges.count(*it) == 0) {
				non_path_edges.push_back(*it);
			}
		}
		return non_path_edges;
	}

	EdgeId FindPathEdge(VertexId v, const set<EdgeId>& path_edges) {
		vector<EdgeId> out_edges = g_.OutgoingEdges(v);
		for (auto it = out_edges.begin(); it != out_edges.end(); ++it) {
			if (path_edges.count(*it) > 0) {
				return *it;
			}
		}
		VERIFY(false);
		return EdgeId(NULL);
	}

	void ProjectComponentPath(const map<VertexId, size_t>& avg_dist,
			const vector<EdgeId>& path, const map<size_t, VertexId>& path_map) {
		DEBUG("Projecting component");
		vector<EdgeId> non_path_edges = NonPathEdges(avg_dist, path);
		set<EdgeId> path_edges_of_all_time(path.begin(), path.end());
		set<size_t> all_dist = AllAvgDist(avg_dist);

		for (auto it = non_path_edges.begin(); it != non_path_edges.end();
				++it) {
			size_t start_dist = avg_dist.find(g_.EdgeStart(*it))->second;
			size_t end_dist = avg_dist.find(g_.EdgeEnd(*it))->second;
			set<size_t> dist_to_split(all_dist.lower_bound(start_dist),
					all_dist.upper_bound(end_dist));
			size_t offset = start_dist;
			EdgeId e = *it;
			for (auto split_it = dist_to_split.begin();
					split_it != dist_to_split.end(); ++split_it) {
				size_t curr = *split_it;
				size_t pos = curr - offset;
				if (pos > 0 && pos < g_.length(e)) {
					DEBUG(
							"Splitting edge " << g_.str(e) << " on position " << pos);
					pair<EdgeId, EdgeId> split_res = g_.SplitEdge(e, pos);
					DEBUG("Splitting edge " << g_.str(e) << " on position " << pos);
					DEBUG("Gluing edges " << g_.str(split_res.first) << " " << g_.str(FindPathEdge(path_map.find(offset)->second,
							path_edges_of_all_time)));
					EdgeId new_edge = g_.GlueEdges(split_res.first,
							FindPathEdge(path_map.find(offset)->second,
									path_edges_of_all_time));
					DEBUG("New edge " << g_.str(new_edge));
					path_edges_of_all_time.insert(new_edge);
					DEBUG(
							"Result: edges " << g_.str(split_res.first) << " " << g_.str(split_res.second));
					e = split_res.second;
					offset = curr;
				}
			}
			path_edges_of_all_time.insert(
					g_.GlueEdges(e,
							FindPathEdge(path_map.find(offset)->second,
									path_edges_of_all_time)));
		}
		DEBUG("Component projected");
	}

	//todo remove
	MappingRange TrivialRange(EdgeId e, size_t& offset) const {
		size_t l = g_.length(e);
		offset += l;
		return MappingRange(Range(offset - l, offset), Range(0, 1));
	}

	MappingPath<EdgeId> TrivialMappingPath(const vector<EdgeId>& edges) const {
		vector<MappingRange> ranges;
		size_t offset = 0;
		for (auto it = edges.begin(); it != edges.end(); ++it) {
			ranges.push_back(TrivialRange(*it, offset));
		}
		return MappingPath<EdgeId>(edges, ranges);
	}

	template<class It>
	void PrintComponent(It begin, It end, size_t cnt) {
		LengthIdGraphLabeler<Graph> labeler(g_);
		WriteComponentsAlongPath(g_, labeler,
				"complex_components/" + ToString(cnt) + ".dot", 5000, 30,
				TrivialMappingPath(vector<EdgeId>(begin, end)),
				*DefaultColorer(g_));
	}
	//end of todo

	void ProcessComponent(const ComponentFinder& comp_finder) {
		static size_t cnt = 0;
		DEBUG("Checking if has conjugate vertices and not single path");
		if (!comp_finder.HaveConjugateVertices()) {
			//find best path!
			pair<size_t, size_t> dist_range = comp_finder.processed().find(
					comp_finder.end_v())->second;

			if (!SinglePath(comp_finder.start_v(), comp_finder.end_v(), dist_range.first,
					dist_range.second)) {
				DEBUG("Check ok");
				DEBUG("Component: " << ++cnt << ". Start vertex - " << g_.int_id(comp_finder.start_v())
						<< ". End vertex - "<<g_.int_id(comp_finder.end_v())<<". Distance ranges "<< dist_range);
				MostCoveredPathChooser<Graph> path_chooser(g_);
				PathProcessor<Graph> best_path_finder(g_, dist_range.first,
						dist_range.second, comp_finder.start_v(), comp_finder.end_v(),
						path_chooser);
				best_path_finder.Process();
				vector<EdgeId> best_path = path_chooser.most_covered_path();
				DEBUG("Best path " << g_.str(best_path));

                remove_dir("complex_components");
                make_dir  ("complex_components");

				PrintComponent(best_path.begin(), best_path.end(), cnt);

				map<VertexId, size_t> dist = AverageDistances(
						comp_finder.processed());

				DEBUG("Licvidating");
				auto split_result = SplitPath(dist, best_path);
				DEBUG("Splitted best path " << split_result.first);
				ProjectComponentPath(dist, split_result.first, split_result.second);
				//						TRACE(
				//								"fake_edge "<<g_.int_id(fake_edge)<< " from "<< g_.int_id(g_.EdgeStart(fake_edge))<<" to "<<g_.int_id(g_.EdgeEnd(fake_edge)));
				Compressor < Graph > (g_).CompressAllVertices();
				//						VertexId v_end = g_.EdgeEnd(fake_edge);
				//						g_.CompressVertex(g_.EdgeStart(fake_edge));
				//						g_.CompressVertex(v_end);

				TRACE("Licvidate finished");
			} else {
				DEBUG("Check fail");
			}
		} else {
			DEBUG("Check fail");
		}
	}

	Graph& g_;
	size_t max_length_;
	size_t length_diff_;

public:
	OppositionLicvidator(Graph& g, size_t max_length, size_t length_diff) :
			g_(g), max_length_(max_length), length_diff_(length_diff) {
	}

	void Licvidate() {
		for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			ComponentFinder comp_finder(g_, max_length_, length_diff_, *it);
			if (comp_finder.TryFindComponent()) {
				ProcessComponent(comp_finder);
			}
		}
	}
private:
	DECL_LOGGER("OppositionLicvidator")
	;
}
;

}
#endif /* BULGE_REMOVER_HPP_ */
