//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <cmath>
#include <stack>
#include <queue>
#include "adt/concurrent_dsu.hpp"
#include "standard_base.hpp"
#include "omni_utils.hpp"
#include "graph_component.hpp"
#include "xmath.h"
#include "sequence/sequence_tools.hpp"
#include "path_processor.hpp"
#include "omni/visualization/visualization.hpp"

namespace omnigraph {

namespace complex_br {

template<class Graph>
class LocalizedComponent: public GraphActionHandler<Graph> /*: public GraphComponent<Graph>*/{
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
	LocalizedComponent(const Graph& g, //It begin, It end,
			VertexId start_vertex/*, const vector<VertexId>& end_vertices*/) :
			base(g, "br_component"), g_(g), start_vertex_(start_vertex) {
		end_vertices_.insert(start_vertex);
		vertex_depth_.insert(make_pair(start_vertex_, Range(0, 0)));
		height_2_vertices_.insert(make_pair(0, start_vertex));
	}

	const Graph& g() const {
		return g_;
	}

	bool IsEndVertex(VertexId v) const {
		FOREACH (EdgeId e, g_.OutgoingEdges(v)) {
			if (contains(g_.EdgeEnd(e)))
				return false;
		}
		return true;
	}

	void AddVertex(VertexId v, Range dist_range) {
//		VERIFY(CheckCloseNeighbour(v));
//		Range r = NeighbourDistanceRange(v);
		DEBUG("Adding vertex " << g_.str(v) << " to the component");
		vertex_depth_.insert(make_pair(v, dist_range));
		height_2_vertices_.insert(make_pair(Average(dist_range), v));
		DEBUG(
				"Range " << dist_range << " Average height " << Average(dist_range));
		FOREACH (EdgeId e, g_.IncomingEdges(v)) {
			end_vertices_.erase(g_.EdgeStart(e));
		}
		if (IsEndVertex(v)) {
			end_vertices_.insert(v);
		}
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
		FOREACH (VertexId v, key_set(vertex_depth_)) {
			if (v == start_vertex_)
				continue;
			if (!AllEdgeIn(v) && !AllEdgeOut(v))
				return false;
		}
		return true;
	}

	bool NeedsProjection() const {
		DEBUG("Checking if component needs projection");
		size_t tot_path_count = TotalPathCount();
		bool answer = tot_path_count > end_vertices_.size();
//		more robust to path processor failure this way VERIFY(tot_path_count >= end_vertices_.size());
		if (answer) {
			DEBUG("Needs projection");
		} else {
			DEBUG("Doesn't need projection");
		}
		return answer;
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
		FOREACH(VertexId v, key_set(vertex_depth_)) {
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
		set<VertexId> vertices = key_set(vertex_depth_);
		return GraphComponent<Graph>(g_, vertices.begin(), vertices.end());
	}

	bool ContainsConjugateVertices() const {
		set<VertexId> conjugate_vertices;
		FOREACH (VertexId v, key_set(vertex_depth_)) {
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
			DEBUG("Deleting vertex " << g_.str(v) << " from the component");
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

	virtual void HandleDelete(EdgeId /*e*/) {
		//empty for now
	}

	virtual void HandleMerge(const vector<EdgeId>& /*old_edges*/, EdgeId /*new_edge*/) {
		VERIFY(false);
	}

	virtual void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) {
		//empty for now
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1, EdgeId /*new_edge_2*/) {
		VertexId start = g_.EdgeStart(old_edge);
		VertexId end = g_.EdgeEnd(old_edge);
		if (contains(start)) {
			VERIFY(vertex_depth_.count(end) > 0);
			VERIFY(avg_distance(end) > avg_distance(start));
			VertexId new_vertex = g_.EdgeEnd(new_edge_1);
			Range new_vertex_depth(distance_range(start));
			new_vertex_depth.shift((int) g_.length(new_edge_1));
			//todo do better later (needs to be synched with splitting strategy)
//					+ (vertex_depth_[end] - vertex_depth_[start])
//							* g_.length(new_edge_1) / g_.length(old_edge);
			DEBUG(
					"Inserting vertex " << g_.str(new_vertex) << " to component during split");
			vertex_depth_.insert(make_pair(new_vertex, new_vertex_depth));
			height_2_vertices_.insert(
					make_pair(Average(new_vertex_depth), new_vertex));
		}
	}

	const multimap<size_t, VertexId>& height_2_vertices() const {
		return height_2_vertices_;
	}

	const set<VertexId> vertices_on_height(size_t height) const {
		set<VertexId> answer;
		for (auto it = height_2_vertices_.lower_bound(height);
				it != height_2_vertices_.upper_bound(height); ++it) {
			answer.insert(it->second);
		}
		return answer;
	}

private:
	DECL_LOGGER("LocalizedComponent")
	;
};

template<class Graph>
class SkeletonTree: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

public:

	const set<EdgeId>& edges() const {
		return edges_;
	}

	const set<VertexId>& vertices() const {
		return vertices_;
	}

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

	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId /*new_edge*/) {
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

	SkeletonTree(const LocalizedComponent<Graph>& br_comp,
			const set<EdgeId>& edges) :
			base(br_comp.g(), "br_tree"), br_comp_(br_comp), edges_(edges) {
		DEBUG("Tree edges " << br_comp.g().str(edges));
		FOREACH(EdgeId e, edges_) {
			vertices_.insert(br_comp_.g().EdgeStart(e));
			vertices_.insert(br_comp_.g().EdgeEnd(e));
		}
	}

private:
	const LocalizedComponent<Graph>& br_comp_;
	set<EdgeId> edges_;
	set<VertexId> vertices_;

private:
	DECL_LOGGER("SkeletonTree")
	;
};

typedef size_t mask;
typedef mask mixed_color_t;
typedef unsigned primitive_color_t;

template<class Graph>
class ComponentColoring: public GraphActionHandler<Graph> {
	typedef GraphActionHandler<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

public:

	size_t CountPrimitiveColors(mixed_color_t color) const {
		size_t cnt = 0;
		for (size_t shift = 0; shift < color_cnt_; ++shift) {
			mixed_color_t prim_color = 1 << shift;
			if ((prim_color & color) != 0) {
				cnt++;
			}
		}
		VERIFY(cnt > 0);
		return cnt;
	}

	primitive_color_t GetAnyPrimitiveColor(mixed_color_t color) const {
		for (size_t shift = 0; shift < color_cnt_; ++shift) {
			if ((1 << shift & color) != 0) {
				return primitive_color_t(shift);
			}
		}
		VERIFY(false);
		return 0;
	}

	bool IsSubset(mixed_color_t super_set, mixed_color_t sub_set) const {
		return (super_set | sub_set) == super_set;
	}

private:

	const LocalizedComponent<Graph>& comp_;
	const size_t color_cnt_;
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
			DEBUG("Coloring exit " << comp_.g().str(v));
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

	ComponentColoring(const LocalizedComponent<Graph>& comp) :
			base(comp.g(), "br_comp_coloring"), comp_(comp), color_cnt_(
					comp_.end_vertices().size()) {
		VERIFY(comp.end_vertices().size() <= sizeof(size_t) * 8);
		ColorComponent();
	}

	mixed_color_t color(VertexId v) const {
		auto it = vertex_colors_.find(v);
		if (it == vertex_colors_.end()) {
			DEBUG("No color for vertex " << comp_.g().str(v));
			DEBUG(
					"Incoming edges " << comp_.g().str(comp_.g().IncomingEdges(v)));
			DEBUG(
					"Outgoing edges " << comp_.g().str(comp_.g().OutgoingEdges(v)));
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

	virtual void HandleMerge(const vector<EdgeId>& /*old_edges*/, EdgeId /*new_edge*/) {
		VERIFY(false);
	}

	virtual void HandleGlue(EdgeId /*new_edge*/, EdgeId edge1, EdgeId edge2) {
		if (comp_.contains(edge1)) {
			VERIFY(comp_.contains(edge2));
			VERIFY(IsSubset(color(edge2), color(edge1)));
		}
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId /*new_edge_2*/) {
		if (comp_.contains(old_edge)) {
			CountAndSetVertexColor(comp_.g().EdgeEnd(new_edge_1));
		}
	}

private:
	DECL_LOGGER("ComponentColoring")
	;
};

template<class Graph>
class SkeletonTreeFinder {

	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef ConcurrentDSU color_partition_ds_t;

	const LocalizedComponent<Graph>& component_;
	const ComponentColoring<Graph>& coloring_;

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
//		VertexId start = component_.g().EdgeStart(e);
		VertexId end = component_.g().EdgeEnd(e);
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
			if (IsGoodEdge(e)) {
				DEBUG("Edge " << component_.g().str(e) << " is classified as good");
				answer.push_back(e);
			} else {
				DEBUG("Edge " << component_.g().str(e) << " is classified as NOT good");
			}
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
		current_level_ = (int) level_heights_.size() - 1;
		size_t end_cnt = 0;
		FOREACH(VertexId v, component_.end_vertices()) {
			good_vertices_.insert(v);
			subtree_coverage_[v] = 0;
			end_cnt++;
		}
	}

	size_t absolute_coverage(EdgeId e) {
		return (size_t) (component_.g().coverage(e) * (double) component_.g().length(e));
	}

	void UpdateNextEdgesAndCoverage(VertexId v) {
		map<mixed_color_t, size_t> best_subtrees_coverage;
		map<mixed_color_t, EdgeId> best_alternatives;
		FOREACH (EdgeId e, component_.g().OutgoingEdges(v)) {
			if (good_edges_.count(e) > 0) {
				VertexId end = component_.g().EdgeEnd(e);
				mixed_color_t color = coloring_.color(e);
				VERIFY(subtree_coverage_.count(end) > 0);
				if (subtree_coverage_[end] + absolute_coverage(e)
						>= best_subtrees_coverage[color]) {
					best_subtrees_coverage[color] = subtree_coverage_[end]
							+ absolute_coverage(e);
					best_alternatives[color] = e;
				}
			}
		}
		size_t coverage = 0;
		FOREACH (size_t cov, value_set(best_subtrees_coverage)) {
			coverage += cov;
		}
		next_edges_[v] = SetAsVector<EdgeId>(value_set(best_alternatives));
		subtree_coverage_[v] = coverage;
	}

public:
    SkeletonTreeFinder(const LocalizedComponent<Graph>& component,
            const ComponentColoring<Graph>& coloring) :
        component_(component),
        coloring_(coloring), 
        level_heights_(SetAsVector<size_t>(component_.avg_distances())), 
        current_level_((int) level_heights_.size() - 1), 
        current_color_partition_(component_.end_vertices().size()) {
        
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
			DEBUG("Processing level " << current_level_ << " on height " << height);
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
						DEBUG("Vertex " << component_.g().str(v) << " is classified as good");
						good_vertices_.insert(v);
						UpdateNextEdgesAndCoverage(v);
					} else {
						DEBUG("Vertex " << component_.g().str(v) << " is classified as NOT good");
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
	DECL_LOGGER("SkeletonTreeFinder")
	;
};

template<class Graph>
void PrintComponent(const LocalizedComponent<Graph>& component,
		const SkeletonTree<Graph>& tree, const string& file_name) {
	typedef typename Graph::EdgeId EdgeId;
	const set<EdgeId> tree_edges = tree.edges();
	shared_ptr<omnigraph::visualization::ElementColorer<typename Graph::EdgeId>> edge_colorer = make_shared<omnigraph::visualization::MapColorer<EdgeId>>(
			tree_edges.begin(), tree_edges.end(),"green", ""
		);
	visualization::WriteComponent(component.AsGraphComponent(), file_name,
			omnigraph::visualization::DefaultColorer(component.g(), edge_colorer),
			*StrGraphLabelerInstance(component.g()));
}

template<class Graph>
void PrintComponent(const LocalizedComponent<Graph>& component,
		const string& file_name) {
	visualization::WriteComponent(component.AsGraphComponent(), file_name,
			omnigraph::visualization::DefaultColorer(component.g()),
			*StrGraphLabelerInstance(component.g()));
}

template<class Graph>
class ComponentProjector {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	const LocalizedComponent<Graph>& component_;
	const ComponentColoring<Graph>& coloring_;
	const SkeletonTree<Graph>& tree_;

//	DEBUG("Result: edges " << g_.str(split_res.first) << " " << g_.str(split_res.second));
//	DEBUG("New vertex" << g_.str(inner_v) << " ");

	void SplitComponent() {
		DEBUG("Splitting component");
		set<size_t> level_heights(component_.avg_distances());
		DEBUG("Level heights " << ToString<size_t>(level_heights));

		GraphComponent<Graph> gc = component_.AsGraphComponent();

		for (auto it = gc.e_begin(); it != gc.e_end(); ++it) {
			VertexId start_v = g_.EdgeStart(*it);
			VertexId end_v = g_.EdgeEnd(*it);
			size_t start_dist = component_.avg_distance(start_v);
			size_t end_dist = component_.avg_distance(end_v);
			DEBUG(
					"Processing edge " << g_.str(*it) << " avg_start " << start_dist << " avg_end " << end_dist);
			set<size_t> dist_to_split(level_heights.lower_bound(start_dist),
					level_heights.upper_bound(end_dist));
			DEBUG("Distances to split " << ToString<size_t>(dist_to_split));

			size_t offset = start_dist;
			EdgeId e = *it;
			for (auto split_it = dist_to_split.begin();
					split_it != dist_to_split.end(); ++split_it) {
				size_t curr = *split_it;
				if (curr == start_dist || curr == end_dist)
					continue;
				DEBUG("Splitting on " << curr);
				size_t pos = curr - offset;
				DEBUG("Splitting edge " << g_.str(e) << " on position " << pos);
				pair<EdgeId, EdgeId> split_res = g_.SplitEdge(e, pos);
				//checks accordance
				VertexId inner_v = g_.EdgeEnd(split_res.first);
				VERIFY(component_.avg_distance(inner_v) == curr);
				e = split_res.second;
				offset = curr;
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
					VERIFY(
							component_.avg_distance(g_.EdgeEnd(e)) == end_height);
					if (tree_.Contains(e)
							&& coloring_.IsSubset(coloring_.color(e), color)) {
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

		for (auto it = SmartSetIterator<Graph, EdgeId>(g_, gc.e_begin(),
				gc.e_end()); !it.IsEnd(); ++it) {
			DEBUG("Trying to project edge " << g_.str(*it));
			EdgeId target = CorrespondingTreeEdge(*it);
			DEBUG("Target found " << g_.str(target));
			if (target != *it) {
				DEBUG(
						"Glueing " << g_.str(*it) << " to target " << g_.str(target));
				g_.GlueEdges(*it, target);
				DEBUG("Glued");
			}
			DEBUG("Edge processed");
		}
		DEBUG("Component projected");
	}

	ComponentProjector(Graph& g, const LocalizedComponent<Graph>& component,
			const ComponentColoring<Graph>& coloring,
			const SkeletonTree<Graph>& tree) :
			g_(g), component_(component), coloring_(coloring), tree_(tree) {

	}

private:
	DECL_LOGGER("ComponentProjector")
	;
};

template<class Graph>
class LocalizedComponentFinder {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	static const size_t exit_bound = 32;
	static const size_t inf = -1ul;

	Graph& g_;
	size_t max_length_;
	size_t length_diff_threshold_;

	LocalizedComponent<Graph> comp_;

	map<VertexId, Range> dominated_;
	set<VertexId> interfering_;

	std::string ToString(EdgeId e) const {
		std::stringstream ss;
		ss << g_.str(e)
				<< " start: "
				<< g_.str(g_.EdgeStart(e))
				<< " end: "
				<< g_.str(g_.EdgeEnd(e));
		return ss.str();
	}

	bool CheckCompleteness() const {
		if (interfering_.size() == 0) {
			VERIFY(comp_.CheckCompleteness());
			return true;
		}
		return false;
	}

	//false if new interfering vertex is not dominated
	//can be slightly modified in new algorithm
	bool ProcessLocality(VertexId processing_v) {
		vector<VertexId> processed_neighb;
		vector<VertexId> unprocessed_neighb;
		FOREACH (EdgeId e, g_.OutgoingEdges(processing_v)) {
			VertexId v = g_.EdgeEnd(e);
			if (!comp_.contains(v)) {
				unprocessed_neighb.push_back(v);
			} else {
				processed_neighb.push_back(v);
			}
		}
		if (!processed_neighb.empty()) {
			FOREACH (VertexId v, unprocessed_neighb) {
				if (dominated_.count(v) > 0) {
					interfering_.insert(v);
				} else {
					return false;
				}
			}
		}
		return true;
	}

	bool AddVertexWithBackwardPaths(VertexId v) {
		DEBUG("Adding vertex with backward paths");
		std::queue<VertexId> q;
		q.push(v);
		while (!q.empty()) {
			VertexId next_v = q.front();
			q.pop();
			if (!ProcessLocality(next_v)) {
				return false;
			}
			if (!comp_.contains(next_v)) {
				VERIFY(dominated_.count(v) > 0);
				comp_.AddVertex(next_v, dominated_.find(next_v)->second);
				FOREACH(EdgeId e, g_.IncomingEdges(next_v)) {
					q.push(g_.EdgeStart(e));
				}
			}
		}
		return true;
	}

	boost::optional<VertexId> ClosestNeigbour() const {
		size_t min_dist = inf;
		boost::optional<VertexId> answer = boost::none;
		for (auto it = dominated_.begin(); it != dominated_.end(); ++it) {
			if (!comp_.contains(it->first) && it->second.start_pos < min_dist) {
				min_dist = it->second.start_pos;
				answer = boost::optional<VertexId>(it->first);
			}
		}
		return answer;
	}

	bool ProcessInterferingVertex(VertexId v) {
		interfering_.erase(v);
		return AddVertexWithBackwardPaths(v);
	}

	bool CheckPathLengths() const {
		VERIFY(CheckCompleteness());
		FOREACH (VertexId v, comp_.end_vertices()) {
			if (comp_.distance_range(v).size() > length_diff_threshold_)
				return false;
		}
		return true;
	}

	bool CheckPositiveHeightDiff() const {
		DEBUG("Checking for positive height diff of each edge");
		GraphComponent<Graph> gc = comp_.AsGraphComponent();
		for (auto it = gc.e_begin(); it != gc.e_end(); ++it) {
			size_t start_height = comp_.avg_distance(g_.EdgeStart(*it));
			size_t end_height = comp_.avg_distance(g_.EdgeEnd(*it));
			//VERIFY(end_height >= start_height);
			if (end_height <= start_height) {
				DEBUG("Check failed for edge " << g_.str(*it) << " start_height " << start_height << " end_height " << end_height);
				return false;
			}
		}
		return true;
	}

	bool CloseComponent() {
		while (!interfering_.empty()) {
			VertexId v = *interfering_.begin();
			DEBUG("Processing interfering vertex " << g_.str(v));
			if (!ProcessInterferingVertex(v)) {
				DEBUG("Vertex processing failed");
				return false;
			}
		}
		return true;
	}

public:
	LocalizedComponentFinder(Graph& g, size_t max_length,
			size_t length_diff_threshold, VertexId start_v) :
			g_(g), max_length_(max_length), length_diff_threshold_(
					length_diff_threshold), comp_(g, start_v) {
		DEBUG(
				"Component finder from vertex " << g_.str(comp_.start_vertex()) << " created");
		DominatedSetFinder<Graph> dominated_set_finder(g_, start_v, max_length);
		dominated_set_finder.FillDominated();
		dominated_ = dominated_set_finder.dominated();
//		ProcessStartVertex();
	}

	bool ProceedFurther() {
		DEBUG("Processing further");

		DEBUG("Choosing closest vertex");
		do {
			optional<VertexId> next_v = ClosestNeigbour();

			if (next_v) {
				DEBUG(
						"Vertex " << g_.str(*next_v) << " was chosen as closest neighbour");
				interfering_.insert(*next_v);
				DEBUG("Trying to construct closure");
				if (!CloseComponent()) {
					DEBUG("Failed to close component");
					return false;
				} else {
					DEBUG("Component closed");
				}
			} else {
				DEBUG("No more vertices can be added");
				return false;
			}
		} while (!comp_.NeedsProjection());

		if (!CheckPathLengths()) {
			DEBUG("Path lengths check failed");
			return false;
		}
		if (!CheckPositiveHeightDiff()) {
			DEBUG("Check for positive height diff of each edge failed");
			return false;
		}
		if (comp_.ContainsConjugateVertices()) {
			DEBUG("Found component contains conjugate vertices");
			return false;
		}
		if (comp_.end_vertices().size() > exit_bound) {
			DEBUG("Too many exits:" << comp_.end_vertices().size());
			return false;
		}
		GraphComponent<Graph> gc = comp_.AsGraphComponent();
		DEBUG("Found component candidate. Vertices: " << g_.str(gc.vertices()));
		return true;
	}

	const LocalizedComponent<Graph>& component() {
		return comp_;
	}

private:
	DECL_LOGGER("LocalizedComponentFinder")
	;
};

template<class Graph>
class ComplexBulgeRemover {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Graph& g_;
	size_t max_length_;
	size_t length_diff_;

	string pics_folder_;

	bool ProcessComponent(LocalizedComponent<Graph>& component,
			size_t candidate_cnt) {
		DEBUG("Processing component");
		ComponentColoring<Graph> coloring(component);
		SkeletonTreeFinder<Graph> tree_finder(component, coloring);
		DEBUG("Looking for a tree");
		if (tree_finder.FindTree()) {
			DEBUG("Tree found");

			SkeletonTree<Graph> tree(component, tree_finder.GetTreeEdges());

			if (!pics_folder_.empty()) {
				PrintComponent(component, tree,
						pics_folder_ + "success/"
								+ ToString(g_.int_id(component.start_vertex()))
								+ "_" + ToString(candidate_cnt) + ".dot");
			}

			ComponentProjector<Graph> projector(g_, component, coloring, tree);
			projector.ProjectComponent();
			DEBUG(
					"Successfully processed component candidate " << candidate_cnt << " start_v " << g_.str(component.start_vertex()));
			return true;
		} else {
			DEBUG(
					"Failed to find skeleton tree for candidate " << candidate_cnt << " start_v " << g_.str(component.start_vertex()));
			if (!pics_folder_.empty()) {
				//todo check if we rewrite all of the previous pics!
				PrintComponent(component,
						pics_folder_ + "fail/"
								+ ToString(g_.int_id(component.start_vertex())) //+ "_" + ToString(candidate_cnt)
								+ ".dot");
			}
			return false;
		}
	}

public:
	ComplexBulgeRemover(Graph& g, size_t max_length, size_t length_diff,
			const string& pics_folder = "") :
			g_(g), max_length_(max_length), length_diff_(length_diff), pics_folder_(
					pics_folder) {
	}

	bool Run() {
        size_t cnt = 0;
		INFO("Complex bulge remover started");
		if (!pics_folder_.empty()) {
//			remove_dir(pics_folder_);
			make_dir(pics_folder_);
			make_dir(pics_folder_ + "success/");
			make_dir(pics_folder_ + "fail/");
		}
		bool something_done_flag = false;
		for (auto it = g_.SmartVertexBegin(); !it.IsEnd(); ++it) {
			DEBUG("Processing vertex " << g_.str(*it));
			size_t candidate_cnt = 0;
			vector<VertexId> vertices_to_post_process;
			{ //important scope!!!
				LocalizedComponentFinder<Graph> comp_finder(g_, max_length_,
						length_diff_, *it);
				while (comp_finder.ProceedFurther()) {
					candidate_cnt++;
					DEBUG(
							"Found component candidate " << candidate_cnt << " start_v " << g_.str(*it));
					LocalizedComponent<Graph> component =
							comp_finder.component();
					if (ProcessComponent(component, candidate_cnt)) {
						something_done_flag = true;
						cnt++;
                        GraphComponent<Graph> gc = component.AsGraphComponent();
						vertices_to_post_process.insert(
								vertices_to_post_process.end(), gc.v_begin(),
								gc.v_end());
						break;
					}
				}
			}
			FOREACH (VertexId v, vertices_to_post_process) {
				it.HandleAdd(v);
				g_.CompressVertex(v);
			}
		}
		INFO("Complex bulge remover finished");
        INFO("Bulges processed " << cnt);
		return something_done_flag;
	}

private:
	DECL_LOGGER("ComplexBulgeRemover")
	;
};

}

}
