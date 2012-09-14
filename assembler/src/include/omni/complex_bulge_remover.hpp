#pragma once

#include <cmath>
#include <stack>
#include <queue>
#include "concurrent_dsu.hpp"
#include "standard_base.hpp"
#include "omni_utils.hpp"
#include "graph_component.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"
#include "path_processor.hpp"

namespace omnigraph {

template<class MapT>
set<typename MapT::key_type> KeySet(const MapT& m) {
	set<typename MapT::key_type> answer;
	for (auto it = m.begin(); it != m.end(); ++it) {
		answer.insert(it->first);
	}
	return answer;
}

template<class MapT>
set<typename MapT::mapped_type> ValueSet(const MapT& m) {
	set<typename MapT::mapped_type> answer;
	for (auto it = m.begin(); it != m.end(); ++it) {
		answer.insert(it->second);
	}
	return answer;
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
					DEBUG(
							"Result: edges " << g_.str(split_res.first) << " " << g_.str(split_res.second) << " inner vertex" << inner_v);
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
					DEBUG(
							"Splitting edge " << g_.str(e) << " on position " << pos);
					DEBUG(
							"Gluing edges " << g_.str(split_res.first) << " " << g_.str(FindPathEdge(path_map.find(offset)->second, path_edges_of_all_time)));
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

			if (!SinglePath(comp_finder.start_v(), comp_finder.end_v(),
					dist_range.first, dist_range.second)) {
				DEBUG("Check ok");
				DEBUG(
						"Component: " << ++cnt << ". Start vertex - " << g_.int_id(comp_finder.start_v()) << ". End vertex - "<<g_.int_id(comp_finder.end_v())<<". Distance ranges "<< dist_range);
				MostCoveredPathChooser<Graph> path_chooser(g_);
				PathProcessor<Graph> best_path_finder(g_, dist_range.first,
						dist_range.second, comp_finder.start_v(),
						comp_finder.end_v(), path_chooser);
				best_path_finder.Process();
				vector<EdgeId> best_path = path_chooser.most_covered_path();
				DEBUG("Best path " << g_.str(best_path));

				remove_dir("complex_components");
				make_dir("complex_components");
				PrintComponent(best_path.begin(), best_path.end(), cnt);

				map<VertexId, size_t> dist = AverageDistances(
						comp_finder.processed());

				DEBUG("Licvidating");
				auto split_result = SplitPath(dist, best_path);
				DEBUG("Splitted best path " << split_result.first);
				ProjectComponentPath(dist, split_result.first,
						split_result.second);
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

//-------------------- new version -------------------

//doesn't support loops
template<class Graph>
class BRComponent: public GraphActionHandler<Graph> /*: public GraphComponent<Graph>*/{

//	typedef GraphComponent<Graph> base;
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	VertexId start_vertex_;
	set<VertexId> end_vertices_;
	//usage of inclusive-inclusive range!!!
	map<VertexId, Range> vertex_depth_;
	multimap<size_t, VertexId> height_2_vertices_;
	size_t diff_threshold_;

	bool AllEdgeOut(VertexId v) const {
		FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
			if (contains(g_.EdgeEnd(e)))
				return false;
		}
		return true;
	}

	bool AllEdgeIn(VertexId v) const {
		FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
			if (!contains(g_.EdgeEnd(e)))
				return false;
		}
		return true;
	}

	size_t Average(Range r) const {
		return r.start_pos;
	}

public:

//	template <class It>
	BRComponent(const Graph& g, //It begin, It end,
			VertexId start_vertex/*, const vector<VertexId>& end_vertices*/) :
			base(g, "br_component"), g_(g), start_vertex_(start_vertex) {
		end_vertices_.insert(start_vertex);
		vertex_depth_.insert(make_pair(start_vertex_, Range(0, 0)));
		height_2_vertices_.insert(make_pair(0, start_vertex));
	}

	const Graph& g() const {
		return g_;
	}

	void AddVertex(VertexId v) {
		DEBUG("Adding vertex " << g_.str(v) << " to the component");
		VERIFY(CheckCloseNeighbour(v));
		Range r = NeighbourDistanceRange(v);
		vertex_depth_.insert(make_pair(v, r));
		height_2_vertices_.insert(make_pair(Average(r), v));
		FOREACH (EdgeId e, g_.IncomingEdges(v)) {
			end_vertices_.erase(g_.EdgeStart(e));
		}
		end_vertices_.insert(v);
	}

	//todo what if path processor will fail inside
	size_t TotalPathCount() const {
		size_t answer = 0;
		FOREACH (VertexId end_v, end_vertices_) {
			PathStorageCallback<Graph> path_storage(g_);
			Range r = vertex_depth_.find(end_v)->second;
			PathProcessor<Graph> best_path_finder(g_, r.start_pos, r.end_pos,
					start_vertex_, end_v, path_storage);
			best_path_finder.Process();
			answer += path_storage.size();
		}
		return answer;
	}

	bool CheckCompleteness() const {
		FOREACH (VertexId v, KeySet(vertex_depth_)) {
			if (!AllEdgeIn(v) && !AllEdgeOut(v))
				return false;
		}
		return true;
	}

	bool NeedsProjection() const {
		size_t tot_path_count = TotalPathCount();
		VERIFY(tot_path_count >= end_vertices_.size());
		return tot_path_count > end_vertices_.size();
	}

	bool contains(VertexId v) const {
		return vertex_depth_.count(v) > 0;
	}

	bool contains(EdgeId e) const {
		return contains(g_.EdgeStart(e)) && contains(g_.EdgeEnd(e));
	}

	Range distance_range(VertexId v) const {
		VERIFY(contains(v));
		return vertex_depth_.find(v)->second;
	}

	size_t avg_distance(VertexId v) const {
		VERIFY(contains(v));
		return Average(vertex_depth_.find(v)->second);
	}

	set<size_t> avg_distances() const {
		set<size_t> distances;
		FOREACH(VertexId v, KeySet(vertex_depth_)) {
			distances.insert(avg_distance(v));
		}
		return distances;
	}

	VertexId start_vertex() const {
		return start_vertex_;
	}

	const set<VertexId>& end_vertices() const {
		return end_vertices_;
	}

	Range NeighbourDistanceRange(VertexId v) const {
		DEBUG("Counting distance range for vertex " << g_.str(v));
		size_t min = numeric_limits<size_t>::max();
		size_t max = 0;
		VERIFY(g_.IncomingEdgeCount(v) > 0);
		VERIFY(CheckCloseNeighbour(v));
		FOREACH (EdgeId e, g_.IncomingEdges(v)) {
			Range range = vertex_depth_.find(g_.EdgeStart(e))->second;
			range.shift(g_.length(e));
			DEBUG("Edge " << g_.str(e) << " provide distance range " << range);
			if (range.start_pos < min)
				min = range.start_pos;
			if (range.end_pos > max)
				max = range.end_pos;
		}
		VERIFY(
				(max > 0) && (min < numeric_limits<size_t>::max()) && (min <= max));
		Range answer(min, max);
		DEBUG("Range " << answer);
		return answer;
	}

	bool CheckCloseNeighbour(VertexId v) const {
		DEBUG("Check if vertex " << g_.str(v) << " can be processed");
		FOREACH (EdgeId e, g_.IncomingEdges(v)) {
			if (!contains(g_.EdgeStart(e))) {
				DEBUG(
						"Blocked by unprocessed or external vertex " << g_.int_id(g_.EdgeStart(e)) << " that starts edge " << g_.int_id(e));
				DEBUG("Check fail");
				return false;
			}
		}
		DEBUG("Check ok");
		return true;
	}

	GraphComponent<Graph> AsGraphComponent() const {
		set<VertexId> vertices = KeySet(vertex_depth_);
		return GraphComponent<Graph>(g_, vertices.begin(), vertices.end());
	}

	bool ContainsConjugateVertices() const {
		set<VertexId> conjugate_vertices;
		FOREACH (VertexId v, KeySet(vertex_depth_)) {
			if (conjugate_vertices.count(v) == 0) {
				conjugate_vertices.insert(g_.conjugate(v));
			} else {
				return true;
			}
		}
		return false;
	}

	virtual void HandleDelete(VertexId v) {
		VERIFY(end_vertices_.count(v) == 0);
		if (contains(v)) {
			DEBUG("Deleting vertex " << g_.str(v) << " from the component=");
			size_t depth = avg_distance(v);
			vertex_depth_.erase(v);
			for (auto it = height_2_vertices_.lower_bound(depth);
					it != height_2_vertices_.upper_bound(depth); ++it) {
				if (it->second == v) {
					height_2_vertices_.erase(it);
					return;
				}
			}
			VERIFY(false);
		}

	}

	virtual void HandleDelete(EdgeId e) {
		//empty for now
	}

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		VERIFY(false);
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		//empty for now
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId new_edge_2) {
		VertexId start = g_.EdgeStart(old_edge);
		VertexId end = g_.EdgeEnd(old_edge);
		if (contains(start)) {
			VERIFY(vertex_depth_.count(end) > 0);
			VERIFY(avg_distance(end) > avg_distance(start));
			VertexId new_vertex = g_.EdgeEnd(new_edge_1);
			Range new_vertex_depth(distance_range(start));
			new_vertex_depth.shift(g_.length(new_edge_1));
			//todo do better later (needs to be synched with splitting strategy)
//					+ (vertex_depth_[end] - vertex_depth_[start])
//							* g_.length(new_edge_1) / g_.length(old_edge);
			DEBUG("Inserting vertex " << g_.str(new_vertex) << " to component during split");
			vertex_depth_.insert(make_pair(new_vertex, new_vertex_depth));
			height_2_vertices_.insert(make_pair(Average(new_vertex_depth), new_vertex));
		}
	}

	const multimap<size_t, VertexId>& height_2_vertices() const {
		return height_2_vertices_;
	}

	const set<VertexId> vertices_on_height(size_t height) const {
		set<VertexId> answer;
		for (auto it = height_2_vertices_.lower_bound(height); it != height_2_vertices_.upper_bound(height); ++it) {
			answer.insert(it->second);
		}
		return answer;
	}

private:
	DECL_LOGGER("BRComponent");
};

template<class Graph>
class BRComponentSpanningTree: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

public:

	bool Contains(EdgeId e) const {
//		VertexId start = br_comp_.g().EdgeStart(e);
//		if (next_edges_.count(start) > 0) {
//			const vector<EdgeId> edges = next_edges_.find(start)->second;
//			return find(e, next_edges_.lower_bound(start), next_edges_.upper_bound(start)) != edges.end();
//		}
//		return false;
		return edges_.count(e) > 0;
	}

	bool Contains(VertexId v) const {
//		return next_edges_.count(v) > 0;
		return vertices_.count(v) > 0;
	}

	virtual void HandleDelete(VertexId v) {
		//verify v not in the tree
		VERIFY(!Contains(v));
	}

	virtual void HandleDelete(EdgeId e) {
		//verify e not in the tree
		DEBUG("Trying to delete " << br_comp_.g().str(e));
		VERIFY(!Contains(e));
	}

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		//verify false
		FOREACH (EdgeId e, old_edges) {
			VERIFY(!Contains(e));
		}
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//		 verify edge2 in tree
//		 put new_edge instead of edge2
		DEBUG("Glueing " << br_comp_.g().str(new_edge) << " " << br_comp_.g().str(edge1) << " " << br_comp_.g().str(edge2));
		if (Contains(edge2)) {
			DEBUG("Erasing from tree: " << br_comp_.g().str(edge2));
			DEBUG("Inserting to tree: " << br_comp_.g().str(new_edge));
			edges_.erase(edge2);
			edges_.insert(new_edge);
		}
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId new_edge_2) {
		if (Contains(old_edge)) {
			edges_.erase(old_edge);
			vertices_.insert(br_comp_.g().EdgeEnd(new_edge_1));
			edges_.insert(new_edge_1);
			edges_.insert(new_edge_2);
		}
	}

	BRComponentSpanningTree(const BRComponent<Graph>& br_comp,
			const set<EdgeId>& edges) :
			base(br_comp.g(), "br_tree"), br_comp_(br_comp), edges_(edges) {
		DEBUG("Tree edges " << br_comp.g().str(edges));
		FOREACH(EdgeId e, edges_) {
			vertices_.insert(br_comp_.g().EdgeStart(e));
			vertices_.insert(br_comp_.g().EdgeEnd(e));
		}
	}

private:
	const BRComponent<Graph>& br_comp_;
	set<EdgeId> edges_;
	set<VertexId> vertices_;

private:
	DECL_LOGGER("BRComponentSpanningTree");
};

typedef size_t mask;
typedef mask mixed_color_t;
typedef unsigned primitive_color_t;

template<class Graph>
class ComponentColoring: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	static const size_t exit_bound = 10;
public:

	static size_t CountPrimitiveColors(mixed_color_t color) {
		size_t cnt = 0;
		for (size_t shift = 0; shift < exit_bound; ++shift) {
			mixed_color_t prim_color = 1 << shift;
			if ((prim_color & color) != 0) {
				cnt++;
			}
		}
		VERIFY(cnt > 0);
		return cnt;
	}

	static primitive_color_t GetAnyPrimitiveColor(mixed_color_t color) {
		for (size_t shift = 0; shift < exit_bound; ++shift) {
			if ((1 << shift & color) != 0) {
				return shift;
			}
		}
		VERIFY(false);
		return 0;
	}

	static bool IsSubset(mixed_color_t super_set, mixed_color_t sub_set) {
		return (super_set | sub_set) == super_set;
	}

private:

	const BRComponent<Graph>& comp_;
	map<VertexId, mixed_color_t> vertex_colors_;

	mixed_color_t CountVertexColor(VertexId v) const {
		mixed_color_t answer = mixed_color_t(0);
		FOREACH(EdgeId e, comp_.g().OutgoingEdges(v)) {
			answer |= color(e);
		}
		return answer;
	}

	void CountAndSetVertexColor(VertexId v) {
		vertex_colors_.insert(make_pair(v, CountVertexColor(v)));
	}

	void ColorComponent() {
		DEBUG("Coloring component");
		size_t cnt = 0;
		FOREACH(VertexId v, comp_.end_vertices()) {
			mixed_color_t color = 1 << cnt;
			vertex_colors_.insert(make_pair(v, color));
			cnt++;
		}
		for (auto it = comp_.height_2_vertices().rbegin();
				it != comp_.height_2_vertices().rend(); ++it) {
			if (vertex_colors_.count(it->second) == 0) {
				DEBUG("Coloring vertex " << comp_.g().str(it->second));
				CountAndSetVertexColor(it->second);
			}
		}
		DEBUG("Component colored");
	}

public:

	ComponentColoring(const BRComponent<Graph>& comp) : base(comp.g(), "br_comp_coloring"), comp_(comp) {
		VERIFY(exit_bound <= sizeof(size_t) * 8);
		VERIFY(comp_.end_vertices().size() < exit_bound);
		ColorComponent();
	}

	mixed_color_t color(VertexId v) const {
		auto it = vertex_colors_.find(v);
		if (it == vertex_colors_.end()) {
			DEBUG("No color for vertex " << comp_.g().str(v));
			DEBUG("Incoming edges " << comp_.g().str(comp_.g().IncomingEdges(v)));
			DEBUG("Outgoing edges " << comp_.g().str(comp_.g().OutgoingEdges(v)));
		}
		VERIFY(it != vertex_colors_.end());
		return it->second;
	}

	mixed_color_t color(EdgeId e) const {
		return color(comp_.g().EdgeEnd(e));
	}

	virtual void HandleDelete(VertexId v) {
		vertex_colors_.erase(v);
	}

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
		VERIFY(false);
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		if (comp_.contains(edge1)) {
			VERIFY(comp_.contains(edge2));
			VERIFY(IsSubset(color(edge2), color(edge1)));
		}
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId new_edge_2) {
		if (comp_.contains(old_edge)) {
			CountAndSetVertexColor(comp_.g().EdgeEnd(new_edge_1));
		}
	}

private:
	DECL_LOGGER("ComponentColoring");
};

template<class Graph>
class BRComponentSpanningTreeFinder {

	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ConcurrentDSU color_partition_ds_t;

	const BRComponent<Graph>& component_;
	const ComponentColoring<Graph>& coloring_;

//	BRComponentSpanningTree<Graph> tree_;
	vector<size_t> level_heights_;

	int current_level_;
	color_partition_ds_t current_color_partition_;

	set<VertexId> good_vertices_;
	set<EdgeId> good_edges_;
	map<VertexId, vector<EdgeId>> next_edges_;
	map<VertexId, size_t> subtree_coverage_;

	bool ConsistentWithPartition(mixed_color_t color) const {
		return current_color_partition_.set_size(
				GetCorrespondingDisjointSet(color))
				== coloring_.CountPrimitiveColors(color);
	}

	bool IsGoodEdge(EdgeId e) const {
		VertexId start = component_.g().EdgeStart(e);
		VertexId end = component_.g().EdgeEnd(e);
		VERIFY(component_.avg_distance(start) < component_.avg_distance(end));
		//check if end is good
		if (good_vertices_.count(end) == 0)
			return false;

//		is subcase of next case
//		//check if end is from previous level
//		if (component_.avg_distance(end) == level_heights_[current_level_+1])
//			return true;

		//check if end color is consistent with partition
		//on level before the start
		return ConsistentWithPartition(coloring_.color(end));
	}

	vector<EdgeId> GoodOutgoingEdges(VertexId v) const {
		vector<EdgeId> answer;
		FOREACH(EdgeId e, component_.g().OutgoingEdges(v)) {
			if (IsGoodEdge(e))
				answer.push_back(e);
		}
		return answer;
	}

	vector<EdgeId> GoodOutgoingEdges(const vector<VertexId>& vertices) const {
		vector<EdgeId> answer;
		FOREACH(VertexId v, vertices) {
			if (component_.end_vertices().count(v) == 0) {
				push_back_all(answer, GoodOutgoingEdges(v));
			}
		}
		return answer;
	}

	set<EdgeId> VectorAsSet(const vector<EdgeId>& edges) const {
		return set<EdgeId>(edges.begin(), edges.end());
	}

	template<class T>
	vector<T> SetAsVector(const set<T>& edges) const {
		return vector<T>(edges.begin(), edges.end());
	}

	primitive_color_t GetCorrespondingDisjointSet(mixed_color_t color) const {
		return current_color_partition_.find_set(
				coloring_.GetAnyPrimitiveColor(color));
	}

	void UpdateColorPartitionWithVertex(VertexId v) {
		VERIFY(component_.g().OutgoingEdgeCount(v) > 0);
		primitive_color_t ds = GetCorrespondingDisjointSet(
				coloring_.color(*(component_.g().OutgoingEdges(v).begin())));
		FOREACH (EdgeId e, component_.g().OutgoingEdges(v)) {
			current_color_partition_.unite(ds,
					GetCorrespondingDisjointSet(coloring_.color(e)));
		}
	}

	bool IsGoodVertex(VertexId v) const {
		if (!ConsistentWithPartition(coloring_.color(v)))
			return false;
		mixed_color_t union_color_of_good_children = mixed_color_t(0);
		FOREACH (EdgeId e, component_.g().OutgoingEdges(v)) {
			if (good_edges_.count(e) > 0) {
				union_color_of_good_children |= coloring_.color(e);
			}
		}
		return coloring_.color(v) == union_color_of_good_children;
	}

	void Init() {
		current_level_ = level_heights_.size() - 1;
		size_t end_cnt = 0;
		FOREACH(VertexId v, component_.end_vertices()) {
			good_vertices_.insert(v);
			subtree_coverage_[v] = 0;
			end_cnt++;
		}
	}

	size_t absolute_coverage(EdgeId e) {
		return component_.g().coverage(e) * component_.g().length(e);
	}

	void UpdateNextEdgesAndCoverage(VertexId v) {
		map<mixed_color_t, size_t> best_subtrees_coverage;
		map<mixed_color_t, EdgeId> best_alternatives;
		FOREACH (EdgeId e, component_.g().OutgoingEdges(v)) {
			if (good_edges_.count(e) > 0) {
				VertexId end = component_.g().EdgeEnd(e);
				mixed_color_t color = coloring_.color(e);
				if (subtree_coverage_[end] + absolute_coverage(e)
						> best_subtrees_coverage[color]) {
					best_subtrees_coverage[color] = subtree_coverage_[end]
							+ absolute_coverage(e);
					best_alternatives[color] = e;
				}
			}
		}
		size_t coverage = 0;
		FOREACH (size_t cov, ValueSet(best_subtrees_coverage)) {
			coverage += cov;
		}
		next_edges_[v] = SetAsVector<EdgeId>(ValueSet(best_alternatives));
		subtree_coverage_[v] = coverage;
	}

public:
	BRComponentSpanningTreeFinder(const BRComponent<Graph>& component,
			const ComponentColoring<Graph>& coloring) :
			component_(component), coloring_(coloring), level_heights_(
					SetAsVector<size_t>(component_.avg_distances())), current_level_(
					level_heights_.size() - 1), current_color_partition_(
					component_.end_vertices().size()) {
		Init();
	}

	const set<EdgeId> GetTreeEdges() const {
		set<EdgeId> answer;
		std::queue<VertexId> vertex_queue;
		vertex_queue.push(component_.start_vertex());
		while (!vertex_queue.empty()) {
			VertexId v = vertex_queue.front();
			vertex_queue.pop();
			if (next_edges_.count(v) == 0)
				continue;
			FOREACH (EdgeId e, next_edges_.find(v)->second) {
				answer.insert(e);
				vertex_queue.push(component_.g().EdgeEnd(e));
			}
		}
		return answer;
	}

	const map<VertexId, vector<EdgeId>>& GetTree() const {
		return next_edges_;
	}

	bool FindTree() {
		DEBUG("Looking for tree");
		while (current_level_ >= 0) {
			size_t height = level_heights_[current_level_];
			set<VertexId> level_vertices = component_.vertices_on_height(
					height);
			VERIFY(!level_vertices.empty());

			//looking for good edges
			insert_all(good_edges_,
					GoodOutgoingEdges(
							vector<VertexId>(level_vertices.begin(),
									level_vertices.end())));

			//counting colors and color partitions
			FOREACH(VertexId v, level_vertices) {
				if (component_.end_vertices().count(v) == 0) {
					UpdateColorPartitionWithVertex(v);
					if (IsGoodVertex(v)) {
						good_vertices_.insert(v);
						UpdateNextEdgesAndCoverage(v);
					}
				}
			}
			current_level_--;
		}
		if (good_vertices_.count(component_.start_vertex()) > 0) {
			DEBUG("Looking for tree was successful");
			return true;
		} else {
			DEBUG("Looking for tree failed");
			return false;
		}
	}

private:
	DECL_LOGGER("BRComponentSpanningTreeFinder");
};

template<class Graph>
class ComponentProjector {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	const BRComponent<Graph>& component_;
	const ComponentColoring<Graph>& coloring_;
	const BRComponentSpanningTree<Graph>& tree_;

	void SplitComponent() {
		DEBUG("Splitting component");
		set<size_t> level_heights(component_.avg_distances());
		GraphComponent<Graph> gc = component_.AsGraphComponent();

		for (auto it = gc.e_begin(); it != gc.e_end(); ++it) {
			VertexId start_v = g_.EdgeStart(*it);
			VertexId end_v = g_.EdgeEnd(*it);
			size_t start_dist = component_.avg_distance(start_v);
			size_t end_dist = component_.avg_distance(end_v);
			set<size_t> dist_to_split(level_heights.lower_bound(start_dist),
					level_heights.upper_bound(end_dist));
			size_t offset = start_dist;
			EdgeId e = *it;
			for (auto split_it = dist_to_split.begin();
					split_it != dist_to_split.end(); ++split_it) {
				size_t curr = *split_it;
				size_t pos = curr - offset;
				if (pos > 0 && pos < g_.length(e)) {
					DEBUG("Splitting edge " << g_.str(e) << " on position " << pos);
					pair<EdgeId, EdgeId> split_res = g_.SplitEdge(e, pos);
					VertexId inner_v = g_.EdgeEnd(split_res.first);
					DEBUG("Result: edges " << g_.str(split_res.first) << " " << g_.str(split_res.second) << " inner vertex" << inner_v);
					//checks accordance
					VERIFY(component_.avg_distance(inner_v) == curr);
					e = split_res.second;
					offset = curr;
				}
			}
		}
		DEBUG("Component split");
	}

	EdgeId CorrespondingTreeEdge(EdgeId e) const {
		DEBUG("Getting height of vertex " << g_.str(g_.EdgeStart(e)));
		size_t start_height = component_.avg_distance(g_.EdgeStart(e));
		DEBUG("Done");
		mixed_color_t color = coloring_.color(e);
		DEBUG("Getting height of vertex " << g_.str(g_.EdgeEnd(e)));
		size_t end_height = component_.avg_distance(g_.EdgeEnd(e));
		DEBUG("Done");
		FOREACH (VertexId v, component_.vertices_on_height(start_height)) {
			if (component_.end_vertices().count(v) == 0) {
				FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
					VERIFY(component_.avg_distance(g_.EdgeEnd(e)) == end_height);
					if (tree_.Contains(e) && coloring_.IsSubset(coloring_.color(e), color)) {
						return e;
					}
				}
			}
		}
		VERIFY(false);
		return EdgeId(NULL);
	}

public:

	void ProjectComponent() {
		SplitComponent();

		DEBUG("Projecting split component");
		GraphComponent<Graph> gc = component_.AsGraphComponent();

		for (auto it = SmartSetIterator<Graph, EdgeId>(g_, gc.e_begin(), gc.e_end()); !it.IsEnd(); ++it) {
			DEBUG("Trying to project edge " << g_.str(*it));
			EdgeId target = CorrespondingTreeEdge(*it);
			DEBUG("Target found " << g_.str(target));
			if (target != *it) {
				DEBUG("Glueing " << g_.str(*it) << " to target " << g_.str(target));
				g_.GlueEdges(*it, target);
				DEBUG("Glued");
			}
			DEBUG("Edge processed");
		}
		DEBUG("Component projected");
	}

	ComponentProjector(Graph& g, const BRComponent<Graph>& component, const ComponentColoring<Graph>& coloring,
			const BRComponentSpanningTree<Graph>& tree) :
			g_(g), component_(component), coloring_(coloring), tree_(tree) {

	}

private:
	DECL_LOGGER("ComponentProjector");
};

template<class Graph>
class BRComponentFinder {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Graph& g_;
	size_t max_length_;
	size_t length_diff_threshold_;

	BRComponent<Graph> comp_;

	set<VertexId> neighbourhood_;
	set<VertexId> can_be_processed_;
	set<VertexId> disturbing_;

	bool CheckCompleteness() const {
		if (disturbing_.size() == 0) {
			VERIFY(comp_.CheckCompleteness());
			return true;
		}
		return false;
	}

	//updating can_be_processed, neighbourhood and disturbing.
	void ProcessLocality(VertexId processing_v) {
		vector<VertexId> processed_neighb;
		vector<VertexId> unprocessed_neighb;
		FOREACH (EdgeId e, g_.OutgoingEdges(processing_v)) {
			VertexId v = g_.EdgeEnd(e);
			if (!comp_.contains(v)) {
				DEBUG("Vertex " << g_.str(v) << " added to neighbourhood")
				neighbourhood_.insert(v);
				if (comp_.CheckCloseNeighbour(v)) {
					DEBUG(
							"Vertex " << g_.str(v) << " added to can_be_processed")
					can_be_processed_.insert(v);
				}
				unprocessed_neighb.push_back(v);
			} else {
				processed_neighb.push_back(v);
			}
		}
		if (!processed_neighb.empty()) {
			FOREACH (VertexId v, unprocessed_neighb) {
				disturbing_.insert(v);
			}
		}
	}

	bool CheckNoEdgeToStart(VertexId v) {
		FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
			if (g_.EdgeEnd(e) == comp_.start_vertex()) {
				return false;
			}
		}
		return true;
	}

	void ProcessStartVertex() {
		DEBUG("Processing start vertex "<< g_.str(comp_.start_vertex()));
		ProcessLocality(comp_.start_vertex());
	}

	void ProcessVertex(VertexId v) {
		//todo delete check
		VERIFY(comp_.CheckCloseNeighbour(v));
		DEBUG("Processing vertex " << g_.str(v));
		comp_.AddVertex(v);
		ProcessLocality(v);

		neighbourhood_.erase(v);
		can_be_processed_.erase(v);
		disturbing_.erase(v);
	}

	bool CheckVertexDist(VertexId v) const {
		return comp_.distance_range(v).start_pos < max_length_;
	}

	bool CheckPathLengths() const {
		VERIFY(CheckCompleteness());
		FOREACH (VertexId v, comp_.end_vertices()) {
			if (comp_.distance_range(v).size() > length_diff_threshold_)
				return false;
		}
		return true;
	}

	VertexId NextVertex() const {
		if (!disturbing_.empty()) {
			return *disturbing_.begin();
		} else {
			VERIFY(!can_be_processed_.empty());
			return *can_be_processed_.begin();
		}
	}

public:
	BRComponentFinder(Graph& g, size_t max_length,
			size_t length_diff_threshold, VertexId start_v) :
			g_(g), max_length_(max_length), length_diff_threshold_(
					length_diff_threshold), comp_(g, start_v) {
	}

	bool TryFindComponent() {
		DEBUG("Trying to find component from vertex " << g_.str(comp_.start_vertex()));
		ProcessStartVertex();
		while (!comp_.CheckCompleteness() || !comp_.NeedsProjection()) {
			if (can_be_processed_.empty()) {
				DEBUG("No more vertices can be processed");
				return false;
			} else {
				VertexId v = NextVertex();
				ProcessVertex(v);
				if (!CheckVertexDist(v) || !CheckNoEdgeToStart(v)) {
					DEBUG(
							"Max component length exceeded or edge to start vertex detected");
					return false;
				}
			}
		}
		if (!CheckPathLengths()) {
			DEBUG("Path lengths check failed");
			return false;
		}
		if (comp_.ContainsConjugateVertices()) {
			DEBUG("Found component contains conjugate vertices");
			return false;
		}
		GraphComponent<Graph> gc = comp_.AsGraphComponent();
		DEBUG("Found component. Vertices: " << g_.str(gc.vertices()));
		return true;
	}

	const BRComponent<Graph>& component() {
		return comp_;
	}

private:
	DECL_LOGGER("BRComponentFinder");
};

template<class Graph>
class ComplexBulgeRemover {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Graph& g_;
	size_t max_length_;
	size_t length_diff_;

	bool ProcessComponent(BRComponent<Graph>& component) {
		ComponentColoring<Graph> coloring(component);
		BRComponentSpanningTreeFinder<Graph> tree_finder(component, coloring);
		if (tree_finder.FindTree()) {
			BRComponentSpanningTree<Graph> tree(component, tree_finder.GetTreeEdges());
			ComponentProjector<Graph> projector(g_, component, coloring, tree);
			projector.ProjectComponent();
			return true;
		} else {
			return false;
		}
	}

public:
	ComplexBulgeRemover(Graph& g, size_t max_length, size_t length_diff) :
			g_(g), max_length_(max_length), length_diff_(length_diff) {
	}

	void Run() {
		for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			vector<VertexId> vertices_to_post_process;
			{
			BRComponentFinder<Graph> comp_finder(g_, max_length_, length_diff_, *it);
			if (comp_finder.TryFindComponent()) {
				BRComponent<Graph> component = comp_finder.component();
				if (ProcessComponent(component)) {
					GraphComponent<Graph> gc = component.AsGraphComponent();
					vertices_to_post_process.insert(vertices_to_post_process.end(),
							gc.v_begin(), gc.v_end());
				}
			}
			}
			FOREACH (VertexId v, vertices_to_post_process) {
				it.HandleAdd(v);
				g_.CompressVertex(v);
			}
		}
	}

private:
	DECL_LOGGER("ComplexBulgeRemover");
};

}
