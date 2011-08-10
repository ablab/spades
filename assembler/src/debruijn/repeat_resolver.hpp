/*
 * repeat_resolver.hpp
 *
 *  Created on: May 5, 2011
 *      Author: antipov
 */

#ifndef REPEAT_RESOLVER_HPP_
#define REPEAT_RESOLVER_HPP_
#include <cmath>
#include <set>
#include <map>
#include <unordered_map>
#include <algorithm>

#include "logging.hpp"
#include "paired_info.hpp"
#include "config.hpp"
#include "omni_utils.hpp"

#include "omni_tools.hpp"
#include "omnigraph.hpp"

#include "ID_track_handler.hpp"
#include "edges_position_handler.hpp"
#include "dijkstra.hpp"

namespace debruijn_graph {

#define MAX_DISTANCE_CORRECTION 10

using omnigraph::SmartVertexIterator;
using omnigraph::Compressor;
using omnigraph::PairedInfoIndex;
using omnigraph::PairInfoIndexData;

template<class Graph>
class RepeatResolver {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	//	typedef SmartVertexIterator<Graph> VertexIter;
	//	typedef omnigraph::SmartEdgeIterator<Graph> EdgeIter;
	typedef omnigraph::PairedInfoIndex<Graph> PIIndex;
	typedef omnigraph::PairInfo<EdgeId> PairInfo;
	typedef vector<PairInfo> PairInfos;

	typedef map<VertexId, set<EdgeId> > NewVertexMap;
	typedef map<VertexId, set<VertexId> > VertexIdMap;

public:

	class EdgeInfo {
	public:
		static const int MAXD = 2;
		static const int MAXSKIPDIST = 4;

		EdgeInfo(const PairInfo &lp_, const int dir_, const EdgeId edge_,
				const int d_) :
				lp(lp_), dir(dir_), edge(edge_), d(d_) {

		}
		inline EdgeId getEdge() {
			return edge;
		}

		bool isClose(int a, int b, double max_diff) {
			return (abs(a - b) < max_diff);
		}
		bool isAdjacent(EdgeInfo other_info, Graph &old_graph,
				Graph &new_graph) {
			//			DEBUG("comparation started: " << edge);
			VertexId v_s = old_graph.EdgeStart(edge);
			VertexId v_e = old_graph.EdgeEnd(edge);
			EdgeId other_edge = other_info.getEdge();
			//			DEBUG("to " << other_edge);
			VertexId other_v_s = old_graph.EdgeStart(other_edge);
			VertexId other_v_e = old_graph.EdgeEnd(other_edge);
			//TODO: insert distance!!!
			int len = old_graph.length(edge);
			int other_len = old_graph.length(other_edge);
			int other_d = other_info.getDistance();

			double max_diff = max(lp.variance, other_info.lp.variance) + 0.5
					+ 1e-9;
			//		max_diff = MAXD;
//ToDo: Understand if it is very dirty hack.
			if ((lp.first != other_info.lp.first)
					&& (new_graph.EdgeStart(lp.first)
							!= new_graph.EdgeEnd(lp.first))
					&& (new_graph.EdgeStart(other_info.lp.first)
							!= new_graph.EdgeEnd(other_info.lp.first))) {
				if ((new_graph.EdgeStart(lp.first)
						== new_graph.EdgeStart(other_info.lp.first))
						|| (new_graph.EdgeEnd(lp.first)
								== new_graph.EdgeEnd(other_info.lp.first)))
					return false;
			}

//TODO:: SHURIK! UBERI ZA SOBOJ !!!
			BoundedDijkstra<Graph, int> dij(old_graph, MAXSKIPDIST);
			dij.run(v_e);
			if (dij.DistanceCounted(other_v_s))
				if (isClose(d + len + dij.GetDistance(other_v_s), other_d,
						max_diff))
					return true;

			dij.run(other_v_e);
			if (dij.DistanceCounted(v_s))
				if (isClose(other_d + other_len + dij.GetDistance(v_s), d,
						max_diff))
					return true;

			if ((other_edge == edge && isClose(d, other_d, max_diff)))
				return true;

			if (lp.first == other_info.lp.first) {
				if ((v_e == other_v_s && isClose(d + len, other_d, max_diff))
						|| (v_s == other_v_e
								&& isClose(d, other_d + other_len, max_diff))
						|| (other_edge == edge && isClose(d, other_d, max_diff))) {
					//				DEBUG("ADJACENT!");
					return true;
				} else {
					//			DEBUG("not adjacent");
					return false;
				}
			}
			return false;
		}

		inline int getDistance() {
			return d;
		}
	public:
		PairInfo lp;
		int dir;
		EdgeId edge;
		int d;

	};

	unordered_map<EdgeId, EdgeId> GetEdgeLabels(){
		return edge_labels;
	}

	RepeatResolver(Graph &old_graph_, IdTrackHandler<Graph> &old_IDs_, int leap,
			PIIndex &ind, EdgesPositionHandler<Graph> &old_pos_,
			Graph &new_graph_, IdTrackHandler<Graph> &new_IDs_,
			EdgesPositionHandler<Graph> &new_pos_) :
			leap_(leap), new_graph(new_graph_), old_graph(old_graph_), new_IDs(
					new_IDs_), old_IDs(old_IDs_), new_pos(new_pos_), old_pos(
					old_pos_) {
		unordered_map<VertexId, VertexId> old_to_new;
		unordered_map<EdgeId, EdgeId> old_to_new_edge;
		cheating_mode = 0;
		cheating_edges.clear();
		size_t paired_size = 0;
		set<VertexId> vertices;
		vertices.clear();
		real_vertices.clear();
		set<EdgeId> edges;
		edges.clear();
		for (auto v_iter = old_graph.SmartVertexBegin(); !v_iter.IsEnd();
				++v_iter) {
			//		if (vertices.find(old_graph.conjugate(*v_iter)) == vertices.end())
			{
				vertices.insert(*v_iter);
				DEBUG(*v_iter);
			}
		}
		for (auto e_iter = old_graph.SmartEdgeBegin(); !e_iter.IsEnd();
				++e_iter) {
			//	if (edges.find(old_graph.conjugate(*e_iter)) == edges.end())
			{
				edges.insert(*e_iter);
			}
		}
		for (auto v_iter = vertices.begin(); v_iter != vertices.end();
				++v_iter) {
			size_t degree = old_graph.IncomingEdgeCount(*v_iter)
					+ old_graph.OutgoingEdgeCount(*v_iter);
			if (degree > 0) {
				VertexId new_vertex = new_graph.AddVertex();
				real_vertices.insert(new_vertex);
				DEBUG("Added vertex" << new_vertex);
				// <<" " << new_graph.conjugate(new_vertex));
				vertex_labels[new_vertex] = *v_iter;
				old_to_new[*v_iter] = new_vertex;
				//	old_to_new[old_graph.conjugate(*v_iter)] = new_graph.conjugate(new_vertex);
			}

		}
		for (auto e_iter = edges.begin(); e_iter != edges.end(); ++e_iter) {
			DEBUG(
					"Adding edge from " << old_to_new[old_graph.EdgeStart(*e_iter)] <<" to " << old_to_new[old_graph.EdgeEnd(*e_iter)]);
			//			DEBUG("Setting coverage to edge length " << old_graph.length(*e_iter) << "  cov: " << old_graph.coverage(*e_iter));
			EdgeId new_edge = new_graph.AddEdge(
					old_to_new[old_graph.EdgeStart(*e_iter)],
					old_to_new[old_graph.EdgeEnd(*e_iter)],
					old_graph.EdgeNucls(*e_iter));
			new_graph.coverage_index().SetCoverage(new_edge,
					old_graph.coverage(*e_iter) * old_graph.length(*e_iter));
//			new_graph.SetCoverage(new_graph.conjugate(new_edge), 0);
			edge_labels[new_edge] = *e_iter;
			DEBUG("Adding edge " << new_edge<< " from" << *e_iter);
			old_to_new_edge[*e_iter] = new_edge;
			//			old_to_new_edge[old_graph.conjugate(*e_iter)] = new_graph.conjugate(new_edge);
			//			PairInfos tmp = ind.GetEdgeInfo(edgeIds[dir][i]);

		}
		old_to_new.clear();
		for (auto p_iter = ind.begin(), p_end_iter = ind.end();
				p_iter != p_end_iter; ++p_iter) {
			PairInfos pi = *p_iter;
			paired_size += pi.size();
			for (size_t j = 0; j < pi.size(); j++) {
				if (old_to_new_edge.find(pi[j].first) != old_to_new_edge.end()
						&& old_to_new_edge.find(pi[j].second)
								!= old_to_new_edge.end()) {
					TRACE(
							"Adding pair " << pi[j].first<<"  " <<old_to_new_edge[pi[j].first] << "  " <<pi[j].second);
					PairInfo *tmp = new PairInfo(old_to_new_edge[pi[j].first],
							pi[j].second, pi[j].d, pi[j].weight, pi[j].variance);
					paired_di_data.AddPairInfo(*tmp, 0);
				} else {
					WARN(
							"Paired Info with deleted edge! " << pi[j].first<<"  " <<pi[j].second);
				}
			}
		}DEBUG("May be size is " << ind.size());
		INFO("paired info size: "<<paired_size);
		assert(leap >= 0 && leap < 100);
	}
	void ResolveRepeats(const string& output_folder);

private:
	int leap_;
	size_t RectangleResolveVertex(VertexId vid);
	size_t CheatingResolveVertex(VertexId vid);

	size_t GenerateVertexPairedInfo(Graph &g, PairInfoIndexData<EdgeId> &ind,
			VertexId vid);
	vector<typename Graph::VertexId> MultiSplit(VertexId v);

	const PairInfo StupidPairInfoCorrector(Graph &new_graph,
			const PairInfo &pair_info) {
		std::multimap<int, EdgeId> Map_queue;
		EdgeId StartEdge = pair_info.first;
		EdgeId EndEdge = pair_info.second;
		int dist = pair_info.d;
		int best = dist + MAX_DISTANCE_CORRECTION + 3;
		//	DEBUG("Adjusting "<<old_IDs.ReturnIntId(edge_labels[StartEdge])<<" "<<old_IDs.ReturnIntId(EndEdge)<<" "<<dist);
		VertexId v;
		vector<EdgeId> edges;
		int len;
		pair<EdgeId, int> Prev_pair;
		if (edge_labels[StartEdge] == EndEdge) {
			if (abs(dist) < MAX_DISTANCE_CORRECTION
			)
				best = 0;
		}
		v = new_graph.EdgeEnd(StartEdge);
		edges = new_graph.OutgoingEdges(v);
		len = new_graph.length(StartEdge);
		for (size_t i = 0; i < edges.size(); i++) {
			//		Prev_pair = make_pair(edges[i], len);
			Map_queue.insert(make_pair(len, edges[i]));
			//		DEBUG("Push ("<<old_IDs.ReturnIntId(edge_labels[edges[i]])<<","<<len<<") ->"<<Map_queue.size());
		}
		while (Map_queue.size() > 0) {
			pair<int, EdgeId> Cur_pair = *(Map_queue.begin());
			//		My_queue.pop();
			Map_queue.erase(Map_queue.begin());
			if (Cur_pair.first - dist < abs(best - dist)) {
				if (edge_labels[Cur_pair.second] == EndEdge) {
					if (abs(Cur_pair.first - dist) < abs(best - dist))
						best = Cur_pair.first;
					//			DEBUG("New best "<<best);
				}
				v = new_graph.EdgeEnd(Cur_pair.second);
				edges.clear();
				edges = new_graph.OutgoingEdges(v);
				len = new_graph.length(Cur_pair.second) + Cur_pair.first;
				for (size_t i = 0; i < edges.size(); i++) {
					//				if ((edges[i] == Prev_pair.second) && (len == Prev_pair.first)) {
					//					DEBUG("SKIP "<<My_queue.size());
					//				} else
					{
						//					Prev_pair = make_pair(edges[i], len);
						typename std::multimap<int, EdgeId>::iterator Map_iter;

						Map_iter = Map_queue.find(len);
						while (Map_iter != Map_queue.end()) {
							if (Map_iter->first != len) {
								Map_iter = Map_queue.end();
								break;
							}
							if (Map_iter->second == edges[i])
								break;
							++Map_iter;
						}

						if (Map_iter == Map_queue.end()) {
							Map_queue.insert(make_pair(len, edges[i]));
							//						DEBUG("Push ("<<edges[i]<<") "<<old_IDs.ReturnIntId(edge_labels[edges[i]])<<","<<len<<") ->"<<Map_queue.size());
						}
					}
				}
			}
		}
		return pair_info.set_distance(best);
	}

	const PairInfo StupidPairInfoCorrectorByOldGraph(Graph &new_graph,
			const PairInfo &pair_info) {
		std::multimap<int, EdgeId> Map_queue;
		EdgeId StartEdge = edge_labels[pair_info.first];
		EdgeId EndEdge = pair_info.second;
		int dist = pair_info.d;
		int best = dist + MAX_DISTANCE_CORRECTION + 3;
		//	DEBUG("Adjusting "<<old_IDs.ReturnIntId(edge_labels[StartEdge])<<" "<<old_IDs.ReturnIntId(EndEdge)<<" "<<dist);
		VertexId v;
		vector<EdgeId> edges;
		int len;
		pair<EdgeId, int> Prev_pair;
		if (StartEdge == EndEdge) {
			if (abs(dist) < MAX_DISTANCE_CORRECTION
			)
				best = 0;
		}
		v = old_graph.EdgeEnd(StartEdge);
		edges = old_graph.OutgoingEdges(v);
		len = old_graph.length(StartEdge);
		for (size_t i = 0; i < edges.size(); i++) {
			//		Prev_pair = make_pair(edges[i], len);
			Map_queue.insert(make_pair(len, edges[i]));
			//		DEBUG("Push ("<<old_IDs.ReturnIntId(edge_labels[edges[i]])<<","<<len<<") ->"<<Map_queue.size());
		}
		while (Map_queue.size() > 0) {
			pair<int, EdgeId> Cur_pair = *(Map_queue.begin());
			//		My_queue.pop();
			Map_queue.erase(Map_queue.begin());
			if (Cur_pair.first - dist < abs(best - dist)) {
				if (Cur_pair.second == EndEdge) {
					if (abs(Cur_pair.first - dist) < abs(best - dist))
						best = Cur_pair.first;
					//			DEBUG("New best "<<best);
				}
				v = old_graph.EdgeEnd(Cur_pair.second);
				edges.clear();
				edges = old_graph.OutgoingEdges(v);
				len = old_graph.length(Cur_pair.second) + Cur_pair.first;
				for (size_t i = 0; i < edges.size(); i++) {
					//				if ((edges[i] == Prev_pair.second) && (len == Prev_pair.first)) {
					//					DEBUG("SKIP "<<My_queue.size());
					//				} else
					{
						//					Prev_pair = make_pair(edges[i], len);
						typename std::multimap<int, EdgeId>::iterator Map_iter;

						Map_iter = Map_queue.find(len);
						while (Map_iter != Map_queue.end()) {
							if (Map_iter->first != len) {
								Map_iter = Map_queue.end();
								break;
							}
							if (Map_iter->second == edges[i])
								break;
							++Map_iter;
						}

						if (Map_iter == Map_queue.end()) {
							Map_queue.insert(make_pair(len, edges[i]));
							//						DEBUG("Push ("<<edges[i]<<") "<<old_IDs.ReturnIntId(edge_labels[edges[i]])<<","<<len<<") ->"<<Map_queue.size());
						}
					}
				}
			}
		}
		if (abs(best - pair_info.d) > 0.000001)
			DEBUG("CORRECTED " << pair_info.d <<" TO " << best);
		PairInfo answer = pair_info;
		answer.d = best;
		return answer;
	}

	pair<bool, PairInfo> CorrectedAndNotFiltered(Graph &new_graph,
			const PairInfo &pair_inf);

	void ResolveEdge(EdgeId eid);
	void dfs(vector<vector<int> > &edge_list, vector<int> &colors, int cur_vert,
			int cur_color);
	VertexIdMap vid_map;
	NewVertexMap new_map;
	Graph &new_graph;
	Graph &old_graph;
	IdTrackHandler<Graph> &new_IDs;
	IdTrackHandler<Graph> &old_IDs;
	EdgesPositionHandler<Graph> &new_pos;
	EdgesPositionHandler<Graph> &old_pos;
	vector<int> edge_info_colors;
	vector<EdgeInfo> edge_infos;
	PairInfoIndexData<EdgeId> paired_di_data;
	unordered_map<VertexId, VertexId> vertex_labels;
	unordered_map<EdgeId, EdgeId> edge_labels;
	set<VertexId> real_vertices;

	int cheating_mode;
	set<EdgeId> cheating_edges;
	map<VertexId, int> resolving_vertices_degrees;
	int sum_count;

private:
	DECL_LOGGER("RepeatResolver")
};

template<class Graph>
vector<typename Graph::VertexId> RepeatResolver<Graph>::MultiSplit(VertexId v) {
	int k = 0;
	for (size_t i = 0; i < edge_info_colors.size(); i++)
		if (edge_info_colors[i] >= k)
			k++;
	DEBUG("splitting to "<< k <<" parts");
	vector<VertexId> res;
	res.resize(k);
	if (k == 1) {
		DEBUG("NOTHING TO SPLIT:( ");
		res[0] = v;
		return res;
	}
	vector<EdgeId> edgeIds[2];
	//TODO: fix labels
	edgeIds[0] = new_graph.OutgoingEdges(v);
	edgeIds[1] = new_graph.IncomingEdges(v);
	vector<unordered_map<EdgeId, EdgeId> > new_edges(k);
	unordered_map<EdgeId, int> old_paired_coverage;
	unordered_map<EdgeId, int> new_paired_coverage;
	//Remember-because of loops there can be same edges in edgeIds[0] and edgeIds[1]
	for (int i = 0; i < k; i++) {
		res[i] = new_graph.AddVertex();
	}DEBUG("Vertex = "<<new_IDs.ReturnIntId(v));
	for (size_t i = 0; i < edge_info_colors.size(); i++) {
		EdgeId le = edge_infos[i].lp.first;
		if (old_paired_coverage.find(le) == old_paired_coverage.end())
			old_paired_coverage[le] = 0;
		old_paired_coverage[le] += edge_infos[i].lp.weight;
		DEBUG(
				"EdgeID = "<<new_IDs.ReturnIntId(le)<<"("<<old_IDs.ReturnIntId(edge_labels[le])<<")"<<" PairEdgeId = "<<old_IDs.ReturnIntId(edge_infos[i].lp.second)<< " Distance = "<<edge_infos[i].lp.d<<" Provided dist = "<<edge_infos[i].getDistance()<<"Weight = "<<edge_infos[i].lp.weight<<" Color = "<<edge_info_colors[i]);
		TRACE(
				"replacing edge " << le<<" with label " << edge_labels[le] << " "<< (new_edges[edge_info_colors[i]].find(le) == new_edges[edge_info_colors[i]].end()));

		EdgeId res_edge = NULL;
		if (edge_infos[i].dir == 0) {
			if (new_graph.EdgeStart(le) != v) {
				ERROR(
						"Non incident edge!!!" << new_graph.EdgeStart(le) <<" instead of "<< v);
			} else {

				if (new_edges[edge_info_colors[i]].find(le)
						== new_edges[edge_info_colors[i]].end())
					res_edge = new_graph.AddEdge(res[edge_info_colors[i]],
							new_graph.EdgeEnd(le), new_graph.EdgeNucls(le));
			}
		} else {
			if (new_graph.EdgeEnd(le) != v) {
				ERROR(
						"Non incident edge!!!" << new_graph.EdgeEnd(le) <<" instead of "<< v);
			} else {
				if (new_edges[edge_info_colors[i]].find(le)
						== new_edges[edge_info_colors[i]].end())
					res_edge = new_graph.AddEdge(new_graph.EdgeStart(le),
							res[edge_info_colors[i]], new_graph.EdgeNucls(le));
			}
		}TRACE("replaced");
		if (res_edge != NULL) {
			new_edges[edge_info_colors[i]].insert(make_pair(le, res_edge));
			edge_labels[res_edge] = edge_labels[le];
			TRACE("before replace first Edge");
			paired_di_data.ReplaceFirstEdge(edge_infos[i].lp, res_edge);

			new_paired_coverage[new_edges[edge_info_colors[i]][le]] = 0;
		} else {
			paired_di_data.ReplaceFirstEdge(edge_infos[i].lp,
					new_edges[edge_info_colors[i]][le]);
		}
		new_paired_coverage[new_edges[edge_info_colors[i]][le]] +=
				edge_infos[i].lp.weight;
		//		old_paired_coverage[le] += edge_infos[i].lp.weight;
	}
	for (int i = 0; i < k; i++) {
		for (auto edge_iter = new_edges[i].begin();
				edge_iter != new_edges[i].end(); edge_iter++) {
			DEBUG(
					"setting coverage to component "<< i << ", from edgeid "<< edge_iter->first<<"("<< new_IDs.ReturnIntId(edge_iter->first)<<")"<<" length "<<new_graph.length(edge_iter->first) <<"   "<< new_graph.coverage(edge_iter->first) <<" taking "<< new_paired_coverage[edge_iter->second] <<"/"<< old_paired_coverage[edge_iter->first]<<" to edgeid "<< edge_iter->second<<"("<< new_IDs.ReturnIntId(edge_iter->second)<<")");
			if ((1.0 * new_paired_coverage[edge_iter->second])
					/ old_paired_coverage[edge_iter->first] > 0.1
					&& (1.0 * new_paired_coverage[edge_iter->second])
							/ old_paired_coverage[edge_iter->first] < 0.9)
				DEBUG("INTERESTING");
			new_graph.coverage_index().SetCoverage(
					edge_iter->second,
					new_graph.length(edge_iter->first)
							* new_graph.coverage(edge_iter->first)
							* new_paired_coverage[edge_iter->second]
							/ old_paired_coverage[edge_iter->first]);
//			new_graph.SetCoverage(new_graph.conjugate(edge_iter->second), 0);
		}
	}

	new_graph.ForceDeleteVertex(v);
	return res;

}

template<class Graph>
void RepeatResolver<Graph>::ResolveRepeats(const string& output_folder) {
	//	old_graph = g;
	//	old_index = ind;
	INFO("resolve_repeats started");
	sum_count = 0;
	for (cheating_mode = 0; cheating_mode < 2; cheating_mode++) {
		bool changed = true;
		set<VertexId> vertices;

		while (changed) {
			changed = false;
			vertices.clear();
			for (auto v_iter = new_graph.SmartVertexBegin(); !v_iter.IsEnd();
					++v_iter) {
				//			if (vertices.find(new_graph.conjugate(*v_iter)) == vertices.end()) {
				vertices.insert(*v_iter);
				//			}
			}
			INFO(
					"Having "<< vertices.size() << "paired vertices, trying to split");
			RealIdGraphLabeler<Graph> IdTrackLabelerAfter(new_graph, new_IDs);
			int GraphCnt = 0;

			omnigraph::WriteSimple(
					output_folder + "resolve_" + ToString(cheating_mode)+"_"+ ToString(GraphCnt) + ".dot",
					"no_repeat_graph", new_graph, IdTrackLabelerAfter);

			for (auto v_iter = vertices.begin(), v_end =
					vertices.end(); v_iter != v_end; ++v_iter) {
				size_t p_size = GenerateVertexPairedInfo(new_graph, paired_di_data, *v_iter);
				DEBUG(" resolving vertex"<<*v_iter<<" "<< p_size);
				int tcount;
				if (cheating_mode == 0)
					tcount = RectangleResolveVertex(*v_iter);
				else
					tcount = CheatingResolveVertex(*v_iter);
				DEBUG("Vertex "<< *v_iter<< " resolved to "<< tcount);
				sum_count += tcount;
				GraphCnt++;
				omnigraph::WriteSimple(
						output_folder + "resolve_" + ToString(cheating_mode)+"_" + ToString(GraphCnt)
								+ ".dot", "no_repeat_graph", new_graph,
						IdTrackLabelerAfter);
			}
		}
	}INFO("total vert" << sum_count);
	INFO("Converting position labels");

	for (auto e_iter = new_graph.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		EdgeId old_edge = edge_labels[*e_iter];
		for (size_t i = 0; i < old_pos.EdgesPositions[old_edge].size(); i++) {
			new_pos.AddEdgePosition(*e_iter,
					old_pos.EdgesPositions[old_edge][i].start_,
					old_pos.EdgesPositions[old_edge][i].end_);
		}
	}
	//	gvis::WriteSimple(  "repeats_resolved_siiimple.dot", "no_repeat_graph", new_graph);
}

template<class Graph>
void RepeatResolver<Graph>::dfs(vector<vector<int> > &edge_list,
		vector<int> &colors, int cur_vert, int cur_color) {
	colors[cur_vert] = cur_color;
	//	DEBUG("dfs-ing, vert num" << cur_vert << " " << cur_color);
	for (int i = 0, sz = edge_list[cur_vert].size(); i < sz; i++) {
		if (colors[edge_list[cur_vert][i]] > -1) {
			if (colors[edge_list[cur_vert][i]] != cur_color) {
				ERROR(
						"error in dfs, neighbour to " << edge_list[cur_vert][i] << " cur_color: "<< cur_color);
			}
		} else {
			if (edge_list[cur_vert][i] != cur_vert) {
				dfs(edge_list, colors, edge_list[cur_vert][i], cur_color);
			}
		}
	}
}

template<class Graph>
pair<bool, PairInfo<typename Graph::EdgeId> > RepeatResolver<Graph>::CorrectedAndNotFiltered(
		Graph &new_graph, const PairInfo &pair_inf) {
	EdgeId right_id = pair_inf.second;
	EdgeId left_id = pair_inf.first;

	if (pair_inf.d - new_graph.length(left_id) > 140) {
		DEBUG(
				"PairInfo "<<edge_labels[left_id]<<"("<<new_graph.length(left_id)<<")"<<" "<<right_id<<"("<<old_graph.length(right_id)<<")"<<" "<<pair_inf.d)DEBUG(
				"too far to correct");
		return make_pair(false, pair_inf);
	}

	PairInfo corrected_info = StupidPairInfoCorrectorByOldGraph(new_graph,
			pair_inf);
	DEBUG(
			"PairInfo "<<edge_labels[left_id]<<" "<<right_id<<" "<<pair_inf.d<< " corrected into "<<corrected_info.d)
	if (abs(corrected_info.d - pair_inf.d) > MAX_DISTANCE_CORRECTION) {
		DEBUG("big correction");
		return make_pair(false, corrected_info);
	}
	if (corrected_info.d - new_graph.length(left_id) > 130) {
		DEBUG("too far");
		return make_pair(false, corrected_info);
	}
	//todo check correctness. right_id belongs to original graph, not to new_graph.
	if (corrected_info.d + new_graph.length(right_id) < 110) {
		DEBUG("too close");
		return make_pair(false, corrected_info);
	}DEBUG("good");
	return make_pair(true, corrected_info);
//	return make_pair(true, pair_inf);
}

template<class Graph>
size_t RepeatResolver<Graph>::GenerateVertexPairedInfo(Graph &new_graph,
		PairInfoIndexData<EdgeId> &paired_data, VertexId vid) {
	DEBUG("Generate vertex paired info for:  " << vid);
	//	DEBUG(new_graph.conjugate(vid));
	edge_infos.clear();
	cheating_edges.clear();
	vector<EdgeId> edgeIds[2];
	edgeIds[0] = new_graph.OutgoingEdges(vid);
	edgeIds[1] = new_graph.IncomingEdges(vid);
	vector<set<EdgeId> > paired_edges;
	DEBUG(edgeIds[0].size()<< "  " << edgeIds[1].size());
	paired_edges.resize(edgeIds[0].size() + edgeIds[1].size());
//TODO:: extract this optional parameter- direction of resolve?
//Or due to r-c structure of original graph it doesn't matter?
	int mult = 1;
	set<EdgeId> right_edges;
	map<VertexId, int> right_vertex_degrees;
	for (int dir = 0; dir < 2; dir++) {
		for (int i = 0, n = edgeIds[dir].size(); i < n; i++) {
			PairInfos tmp = paired_di_data.GetEdgeInfos(edgeIds[dir][i]);
			for (int j = 0, sz = tmp.size(); j < sz; j++) {
				EdgeId right_id = tmp[j].second;
				EdgeId left_id = tmp[j].first;
				double d = tmp[j].d;
				//				int w = tmp[j].weight;
				//				if (w < 10) continue;
				int dif_d = 0;
				//				if ((d >=new_graph.length(left_id))||(edge_labels[left_id] == right_id ))
				{
					if ((dir == 1) /*&& (edge_labels[left_id] != right_id)*/) {
						dif_d = new_graph.length(left_id);

					}
					if (d * mult >= -0.001) {
//						if (cheating_mode && (i == 1) && (right_id == left_id)
//								&& (tmp[j].d == 0)) {
//							DEBUG("Paired info in cheating mode ignored");
//							//ignoring information from incoming edge to itself, ignoring
//							cheating_edges.insert(left_id);
//							continue;
//						}

						DEBUG("PairInfo: " << new_IDs.ReturnIntId(tmp[j].first)<<" "<<old_IDs.ReturnIntId(edge_labels[tmp[j].first]) << " " << old_IDs.ReturnIntId(tmp[j].second) <<" "<< tmp[j].d);

						pair<bool, PairInfo> correction_result =
								CorrectedAndNotFiltered(new_graph, tmp[j]);
						if (!correction_result.first)
							continue;
						DEBUG(
								"PairInfo "<<edge_labels[left_id]<<" "<<right_id<<" "<<d<< " corrected into "<<tmp[j].d<< "weight" << tmp[j].weight);
						DEBUG(
								"PairInfo: " << old_IDs.ReturnIntId(edge_labels[tmp[j].first]) << " " << old_IDs.ReturnIntId(tmp[j].second) <<" "<< tmp[j].d);
						EdgeInfo ei(correction_result.second, dir, right_id,
								correction_result.second.d - dif_d);
						edge_infos.push_back(ei);
						//					DEBUG(right_id);
						right_edges.insert(right_id);
					}

				}
			}
		}
	}
	return right_edges.size();
}

template<class Graph>
size_t RepeatResolver<Graph>::RectangleResolveVertex(VertexId vid) {
	DEBUG("Rectangle resolve vertex started");
	int size = edge_infos.size();
	edge_info_colors.resize(size);
	for (int i = 0; i < size; i++) {
		edge_info_colors[i] = -1;
	}
	vector<vector<int> > neighbours;
	neighbours.resize(size);
	cheating_edges.clear();

	DEBUG("constructing edge-set");
	set<EdgeId> edges;
	edges.clear();

	for (auto e_iter = old_graph.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		{
			edges.insert(*e_iter);
		}
	}DEBUG("checking edge-infos ");
	for (int i = 0; i < size; i++) {
		if (edges.find(edge_infos[i].getEdge()) == edges.end()) {
			ERROR("fake edge: "<< edge_infos[i].getEdge());
		}
	}
	for (int i = 0; i < size; i++)
		neighbours[i].resize(0);
	for (int i = 0; i < size; i++) {
		//		DEBUG("Filling " << i << " element");
		if (edges.find(edge_infos[i].getEdge()) == edges.end()) {
			ERROR("fake edge");
		}
		for (int j = 0; j < size; j++) {
			if (edge_infos[i].isAdjacent(edge_infos[j], old_graph, new_graph)
					&& !edge_infos[j].isAdjacent(edge_infos[i], old_graph,
							new_graph))
				WARN(
						"ASSYMETRIC: " << new_IDs.ReturnIntId(edge_infos[i].getEdge()) << " " << new_IDs.ReturnIntId(edge_infos[j].getEdge()));
			if (edge_infos[i].isAdjacent(edge_infos[j], old_graph, new_graph)) {
				neighbours[i].push_back(j);
				neighbours[j].push_back(i);
				DEBUG(
						old_IDs.ReturnIntId(edge_infos[i].lp.second) <<" " << edge_infos[i].d << "is adjacent "<<old_IDs.ReturnIntId( edge_infos[j].lp.second) <<" " << edge_infos[j].d);

			}
		}
	}
	int cur_color = 0;
	DEBUG("dfs started");
	for (int i = 0; i < size; i++) {
		if (edge_info_colors[i] == -1) {
			dfs(neighbours, edge_info_colors, i, cur_color);
			cur_color++;
		}
	}
	MultiSplit(vid);
	return cur_color;
}


template<class Graph>
size_t RepeatResolver<Graph>::CheatingResolveVertex(VertexId vid) {
	DEBUG("ACHTUNG, cheating resolve vertex started "<<new_IDs.ReturnIntId(vid));
	int size = edge_infos.size();
	edge_info_colors.resize(size);
	vector<EdgeId> edgeIds[2];
	edgeIds[0] = new_graph.OutgoingEdges(vid);
	edgeIds[1] = new_graph.IncomingEdges(vid);
	map<EdgeId, int> EdgeIdMap[2];
	size_t out_count = edgeIds[0].size();
	size_t in_count = edgeIds[1].size();
	size_t counts[2];
	counts[0] = edgeIds[0].size();
	counts[1] = edgeIds[1].size();
	for(int ind = 0; ind < 2; ind ++)
		for(size_t i = 0; i < counts[ind]; i++) {
			DEBUG("direction" << ind << " edge " << new_IDs.ReturnIntId(edgeIds[ind][i]));
			EdgeIdMap[ind].insert(make_pair(edgeIds[ind][i], i));
		}
//	int table[out_count][in_count];
//	for(int i = 0; i < out_count; i ++)
//		for(int j = 0; j < in_count; j++)
//			table[i][j] = 0;
	vector<vector<int> > neighbours;
	neighbours.resize(in_count + out_count);
	for(int i = 0; i < size; i++){
		DEBUG("info N "<<i<<":" <<new_IDs.ReturnIntId(edge_infos[i].lp.first)<<" -> "<<old_IDs.ReturnIntId(edge_infos[i].lp.second)<<" dist "<< edge_infos[i].d);
	}
	for(int i = 0; i < size; i++){
		EdgeId second = NULL;
		EdgeId first = edge_infos[i].lp.first;
		DEBUG("trying first "<< new_IDs.ReturnIntId(first)<<" with paired "<< old_IDs.ReturnIntId(edge_infos[i].lp.second));
		if (EdgeIdMap[0].find(first) == EdgeIdMap[0].end())
			continue;
		for(int j = 0; j < size; j++){
			if ( (first != edge_infos[j].lp.first)&&(edge_infos[i].d==edge_infos[j].d)/*(edge_infos[i].isAdjacent(edge_infos[j], old_graph, new_graph))*/
				&& (edge_infos[i].lp.second == edge_infos[j].lp.second)	){
				if (second == NULL) {
					second = edge_infos[j].lp.first;
				}
				else {
					if (second != edge_infos[j].lp.first){
						second = NULL;
						DEBUG("multiple pairing, break");
						break;
					}
				}

			}
		}
		if (second != NULL) {
			DEBUG("found second "<< new_IDs.ReturnIntId(second));

			if (EdgeIdMap[1].find(second) == EdgeIdMap[1].end())
				continue;
			else {
				size_t first_ind = EdgeIdMap[0][first];
				size_t second_ind = EdgeIdMap[1][second];
				DEBUG(first_ind <<" "<< second_ind);
				neighbours[first_ind + counts[1]].push_back(second_ind);
				neighbours[second_ind].push_back(first_ind + counts[1]);
				DEBUG("neighbours "<<first_ind<<" + "<<counts[0]<<"  "<<second_ind);
//				table[first_ind][second_ind] = 1;
			}
		}
	}
	DEBUG("cheater_colors creating");
	vector<int> cheater_colors(counts[0] + counts[1], -1);

	int cur_color = 0;
	DEBUG("dfs started");
	for (size_t i = 0; i < counts[0] + counts[1]; i++) {
		if (cheater_colors[i] == -1) {
			dfs(neighbours, cheater_colors, i, cur_color);
			cur_color++;
		}
	}
	DEBUG("Colours "<<cheater_colors);

	bool bad = true;
	for (size_t i = 0; i < counts[0] + counts[1]; i++) {
		bad = true;
		for (size_t j = 0; j < counts[0] + counts[1]; j++) {
			if ((i != j) && (cheater_colors[i] == cheater_colors[j])){
				bad = false;
				break;
			}
		}
		if (bad) break;
	}

	if (bad) {
		DEBUG("Cheat failed");
		return 1;
	}

	for (size_t i = 0; i<edge_info_colors.size();i++){
		if (EdgeIdMap[0].find(edge_infos[i].lp.first) != EdgeIdMap[0].end()){
			edge_info_colors[i] = cheater_colors[EdgeIdMap[0][edge_infos[i].lp.first]+ counts[1]];
		} else {
			if (EdgeIdMap[1].find(edge_infos[i].lp.first) != EdgeIdMap[1].end()){
				edge_info_colors[i] = cheater_colors[EdgeIdMap[1][edge_infos[i].lp.first] ];
			}
		}
	}

	MultiSplit(vid);
	return cur_color;
}

}

#endif /* REPEAT_RESOLVER_HPP_ */
