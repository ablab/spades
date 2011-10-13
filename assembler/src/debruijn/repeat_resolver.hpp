/*
 * Repeat_resolver.hpp
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
#include "simple_tools.hpp"
#include "omni/paired_info.hpp"
#include "deleted_vertex_handler.hpp"
#include "config_struct.hpp"
#include "omni/omni_utils.hpp"

#include "omni/omni_tools.hpp"
#include "omni/omnigraph.hpp"
#include "omni/splitters.hpp"

#include "omni/ID_track_handler.hpp"
#include "omni/edges_position_handler.hpp"
#include "omni/total_labeler.hpp"
#include "omni/dijkstra.hpp"
#include "new_debruijn.hpp"


namespace debruijn_graph {

#define MAX_DISTANCE_CORRECTION 10



using omnigraph::SmartVertexIterator;
using omnigraph::Compressor;
using omnigraph::PairedInfoIndex;
using omnigraph::PairInfoIndexData;
using debruijn_graph::DeletedVertexHandler;
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



		bool IsEdgesOnDistanceAdjacent(EdgeId edge,       int d,
									   EdgeId other_edge, int other_d,
									   Graph &old_graph, double max_diff, bool first_equal ){

			VertexId v_s = old_graph.EdgeStart(edge);
			VertexId v_e = old_graph.EdgeEnd(edge);
//			EdgeId other_edge = other_info.getEdge();
			//			DEBUG("to " << other_edge);
			VertexId other_v_s = old_graph.EdgeStart(other_edge);
			VertexId other_v_e = old_graph.EdgeEnd(other_edge);
			//TODO: insert distance!!!
			int len = old_graph.length(edge);
			int other_len = old_graph.length(other_edge);
//			int other_d = other_info.getDistance();




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

			if (first_equal) {
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





		bool isAdjacent(EdgeInfo other_info, Graph &old_graph,
				Graph &new_graph, EdgeLabelHandler<Graph> &labels_after) {
			//			DEBUG("comparation started: " << edge);
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
			double max_diff = max(lp.variance, other_info.lp.variance) + 0.5
					+ 1e-9;

			bool old_res = IsEdgesOnDistanceAdjacent(this->edge, this->d, other_info.getEdge()
					,other_info.getDistance(), old_graph, max_diff, lp.first == other_info.lp.first);


			//new_version
			set<EdgeId> edges_set = labels_after.edge_inclusions[this->edge];
			set<EdgeId> other_edges_set = labels_after.edge_inclusions[other_info.getEdge()];

			bool new_res = false;


			for(auto this_edge_it = edges_set.begin(); this_edge_it != edges_set.end(); ++ this_edge_it)
				for(auto other_edge_it = other_edges_set.begin(); other_edge_it != other_edges_set.end(); ++ other_edge_it)
					if( IsEdgesOnDistanceAdjacent(*this_edge_it, this->d, *other_edge_it
							,other_info.getDistance(), new_graph, max_diff, lp.first == other_info.lp.first))
					new_res  = true;

			if (old_res != new_res) {
				DEBUG("difference in isAdjacent for ("<<this->getEdge()<<", ("<<this->lp.first<<", "<<this->lp.second<<", "<<this->lp.d<<"), "<<this->d<<")");
				DEBUG("                          VS ("<<other_info.getEdge()<<", ("<<other_info.lp.first<<", "<<other_info.lp.second<<", "<<other_info.lp.d<<"), "<<other_info.d<<")");
				DEBUG("   old is "<<old_res<<"    new is "<<new_res);
				DEBUG("   first set size "<<edges_set.size()<<"    second set size "<<other_edges_set.size());
			}
			return old_res;
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
			EdgesPositionHandler<Graph> &new_pos_,
			DeletedVertexHandler<Graph> &deleted_handler_,
			EdgeLabelHandler<Graph> &LabelsAfter_
			) :
			leap_(leap), new_graph(new_graph_), old_graph(old_graph_), new_IDs(
					new_IDs_), old_IDs(old_IDs_), new_pos(new_pos_), old_pos(
					old_pos_), deleted_handler(deleted_handler_),  labels_after(LabelsAfter_){
		TRACE("Constructor started");
		unordered_map<VertexId, VertexId> old_to_new;
		unordered_map<EdgeId, EdgeId> old_to_new_edge;
		cheating_mode = 0;
		rc_mode = cfg::get().rr.symmetric_resolve;
		global_cheating_edges.clear();
		size_t paired_size = 0;
		set<VertexId> vertices;
		set<VertexId> rc_vertices;
		vertices.clear();
		real_vertices.clear();
		set<EdgeId> edges;
		set<EdgeId> rc_edges;
		edges.clear();
		near_vertex = cfg::get().rr.near_vertex;
		for (auto v_iter = old_graph.SmartVertexBegin(); !v_iter.IsEnd();
				++v_iter) {
			{
				vertices.insert(*v_iter);
				TRACE(*v_iter);
			}
		}
		for (auto e_iter = old_graph.SmartEdgeBegin(); !e_iter.IsEnd();
				++e_iter) {
			{
				edges.insert(*e_iter);
				TRACE("edge added to array " << *e_iter);
			}
		}
		for (auto v_iter = vertices.begin(); v_iter != vertices.end();
				++v_iter) {
			if (rc_mode) {
				if (rc_vertices.find(*v_iter) == rc_vertices.end())
					rc_vertices.insert(conj_wrap(old_graph, *v_iter));
				else
					continue;
			}

			size_t degree = old_graph.IncomingEdgeCount(*v_iter)
					+ old_graph.OutgoingEdgeCount(*v_iter);
			if (degree > 0) {
				VertexId new_vertex = new_graph.AddVertex();
//				new_IDs.AddVertexIntId(new_vertex, old_IDs.ReturnIntId(*v_iter));
				real_vertices.insert(new_vertex);

				TRACE("Added vertex" << new_vertex);
				vertex_labels[new_vertex] = *v_iter;
				old_to_new[*v_iter] = new_vertex;
				if (rc_mode) {
					VertexId new_rc_vertex = conj_wrap(new_graph, new_vertex);
					VertexId old_rc_vertex = conj_wrap(old_graph, *v_iter);
					real_vertices.insert(new_rc_vertex);
					vertex_labels[new_rc_vertex] = old_rc_vertex;
					old_to_new[old_rc_vertex] = new_rc_vertex;
				}
			}

		}
		for (auto e_iter = edges.begin(); e_iter != edges.end(); ++e_iter) {
			if (rc_mode) {
				if (rc_edges.find(*e_iter) == rc_edges.end())
					rc_edges.insert(conj_wrap(old_graph, *e_iter));
				else
					continue;
			}
			TRACE(
					"Adding edge from " << old_to_new[old_graph.EdgeStart(*e_iter)] <<" to " << old_to_new[old_graph.EdgeEnd(*e_iter)]);
			//			DEBUG("Setting coverage to edge length " << old_graph.length(*e_iter) << "  cove: " << old_graph.coverage(*e_iter));
			EdgeId new_edge = new_graph.AddEdge(
					old_to_new[old_graph.EdgeStart(*e_iter)],
					old_to_new[old_graph.EdgeEnd(*e_iter)],
					old_graph.EdgeNucls(*e_iter));
//			new_IDs.AddEdgeIntId(new_edge, old_IDs.ReturnIntId(*e_iter));
			WrappedSetCoverage(new_edge,
					old_graph.coverage(*e_iter) * old_graph.length(*e_iter));

			edge_labels[new_edge] = *e_iter;
			TRACE("Adding edge " << new_edge<< " from" << *e_iter);
			old_to_new_edge[*e_iter] = new_edge;
			new_pos.AddEdgePosition(new_edge, old_pos.edges_positions().find(*e_iter)->second);


			if (rc_mode) {
				EdgeId new_rc_edge = conj_wrap(new_graph, new_edge);
				EdgeId old_rc_edge = conj_wrap(old_graph, *e_iter);
				edge_labels[new_rc_edge] = old_rc_edge;
				old_to_new_edge[old_rc_edge] = new_rc_edge;
				new_pos.AddEdgePosition(new_rc_edge, old_pos.edges_positions().find(old_rc_edge)->second);
				TRACE("rc edge added");
			}

		}
		TRACE("Edge Adding finished");
		old_to_new.clear();
		DEBUG("Copying of paired info started");
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
		VERIFY(leap >= 0 && leap < 100);
	}
	void ResolveRepeats(const string& output_folder);

private:
	int leap_;
	int near_vertex;
	map<int, VertexId> fillVerticesAuto();
	map<int, VertexId> fillVerticesComponents();
	size_t RectangleResolveVertex(VertexId vid);
	size_t CheatingResolveVertex(VertexId vid);
	void BanRCVertex(VertexId v );
	ConjugateDeBruijnGraph::VertexId conj_wrap(ConjugateDeBruijnGraph& g, ConjugateDeBruijnGraph::VertexId v);
	NonconjugateDeBruijnGraph::VertexId conj_wrap(NonconjugateDeBruijnGraph& g, NonconjugateDeBruijnGraph::VertexId v);

	ConjugateDeBruijnGraph::EdgeId conj_wrap(ConjugateDeBruijnGraph& g, ConjugateDeBruijnGraph::EdgeId e);
	NonconjugateDeBruijnGraph::EdgeId conj_wrap(NonconjugateDeBruijnGraph& g, NonconjugateDeBruijnGraph::EdgeId e);
	bool rc_mode;
	void WrappedSetCoverage(EdgeId e, int cov);
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
	DeletedVertexHandler<Graph> &deleted_handler;
	EdgeLabelHandler<Graph> &labels_after;
	vector<int> edge_info_colors;
	vector<EdgeInfo> edge_infos;
	PairInfoIndexData<EdgeId> paired_di_data;
	unordered_map<VertexId, VertexId> vertex_labels;
	unordered_map<EdgeId, EdgeId> edge_labels;
	set<VertexId> real_vertices;


	int cheating_mode;
	map<EdgeId, int> local_cheating_edges;
	set<EdgeId> global_cheating_edges;
	map<VertexId, int> resolving_vertices_degrees;
	int sum_count;

private:
	DECL_LOGGER("RepeatResolver")
};

template<class Graph>
vector<typename Graph::VertexId> RepeatResolver<Graph>::MultiSplit(VertexId v) {
	int k = 0;
	vector<EdgeId> edgeIds[2];
	//TODO: fix labels
	edgeIds[0] = new_graph.OutgoingEdges(v);
	edgeIds[1] = new_graph.IncomingEdges(v);
	map<EdgeId, int> edgeCounts;
	for(int i = 0; i < 2; i++) {
		for(size_t j = 0; j < edgeIds[i].size(); j++)
			edgeCounts.insert(make_pair(edgeIds[i][j], 0));
	}

	for (size_t i = 0; i < edge_info_colors.size(); i++) {
		if (edge_info_colors[i] >= k)
			k = edge_info_colors[i];
		EdgeId le = edge_infos[i].lp.first;
		edgeCounts[le] ++;
	}
	vector<VertexId> res;
	if (k == 0) {
		DEBUG("NOTHING TO SPLIT:( ");
//		for (size_t j = 0; j < edge_infos.size(); j++){
//			if (edge_info_colors[j] == 1)
//				paired_di_data.ReplaceFirstEdge(edge_infos[j].lp, edge_infos[j].lp.first);
//		}
		res.resize(1);
		res[0] = v;
		return res;
	}

	for(auto iter = edgeCounts.begin(); iter != edgeCounts.end(); ++iter) {
		if (iter->second == 0 && cheating_mode == 2) {
			INFO("Adding no-paired edge: " << new_IDs.ReturnIntId(iter->first)<< " potential bug here.");
			PairInfos tmp = paired_di_data.GetEdgeInfos(iter->first);
			int fl = 0;
			for(size_t j = 0; j < tmp.size(); j ++ ){
				EdgeId right_id = tmp[j].second;
//				EdgeId left_id = tmp[j].first;
				double d = tmp[j].d;
				int w = tmp[j].weight;
				if (w < 1e-8)
					continue;
			//	if (v == new_graph.S)
				int dif_d = 0;
				fl ++;
				int dir = 0;
				// it doesn't matter now
				EdgeInfo ei(tmp[j], dir, right_id,
					int(d - dif_d));
				edge_infos.push_back(ei);
				edge_info_colors.push_back(k);
			}
			if (fl)
			k++;
		}
	}
	k++;
	DEBUG("splitting to "<< k <<" parts");

	res.resize(k);
	if (k == 1) {
		DEBUG("NOTHING TO SPLIT:( ");
//		for (size_t j = 0; j < edge_infos.size(); j++){
//			if (edge_info_colors[j] == 1)
//				paired_di_data.ReplaceFirstEdge(edge_infos[j].lp, edge_infos[j].lp.first);
//		}
		res[0] = v;
		return res;
	}


	for(auto iter = edgeCounts.begin(); iter != edgeCounts.end(); ++iter) {
		if (iter->second > 1) {
			paired_di_data.DeleteEdgeInfo(iter->first);
		} else 
		if (iter->second == 1){
			int updated_edge_color = -1;
			for (size_t j = 0; j < edge_infos.size(); j++){
				if (edge_infos[j].lp.first == iter->first){
					if (updated_edge_color == -1){
						updated_edge_color = edge_info_colors[j];
					} else {
						if (updated_edge_color != edge_info_colors[j]){
							WARN("Different colors found for one colored edge info");
						}
					}
				}
			}

			if (updated_edge_color > -1){
				for (size_t j = 0; j < edge_infos.size(); j++){
					if ((edge_info_colors[j] == updated_edge_color)&&(edge_infos[j].lp.first == iter->first)){
						edge_info_colors.erase(edge_info_colors.begin()+j);
						edge_infos.erase(edge_infos.begin()+j);
						j--;
					}
				}
				PairInfos tmp = paired_di_data.GetEdgeInfos(iter->first);
				for(size_t info_j = 0; info_j < tmp.size(); info_j ++ ){
					EdgeInfo ei(tmp[info_j], 0, tmp[info_j].second, 0); //check
					edge_infos.push_back(ei);
					edge_info_colors.push_back(updated_edge_color);
				}
				paired_di_data.DeleteEdgeInfo(iter->first);
			}
		}
	}
	
	vector<unordered_map<EdgeId, EdgeId> > new_edges(k);
	vector<vector<EdgeId> > edges_for_split(k);

	unordered_map<EdgeId, int> old_paired_coverage;
	vector<unordered_map<EdgeId, int>> colored_paired_coverage(k);
	//Remember-because of loops there can be same edges in edgeIds[0] and edgeIds[1]
//TODO: fix loop problem by split loop into 2 eges.

	//	for (int i = 0; i < k; i++) {
//		res[i] = new_graph.AddVertex();
//	}DEBUG("Vertex = "<<new_IDs.ReturnIntId(v));

	for (size_t i = 0; i < edge_info_colors.size(); i++) {
		EdgeId le = edge_infos[i].lp.first;
		int color = edge_info_colors[i];
		if (old_paired_coverage.find(le) == old_paired_coverage.end())
			old_paired_coverage[le] = 0;
		old_paired_coverage[le] += edge_infos[i].lp.weight;
		if (colored_paired_coverage[color].find(le) == colored_paired_coverage[color].end())
			colored_paired_coverage[color][le] = 0;
		colored_paired_coverage[color][le] += edge_infos[i].lp.weight;
	}

	for (int i = 0; i < k; i++) {
		vector<EdgeId> split_edge;
		vector<double> split_coeff;
		for (auto it = colored_paired_coverage[i].begin(); it !=  colored_paired_coverage[i].end(); it++){
			if (it->second != 0){
//		if (((double)it->second)/old_paired_coverage[it->first] > 0.01){
					split_edge.push_back(it->first);
					if (local_cheating_edges.find(it->first) != local_cheating_edges.end()) {
						DEBUG("local_cheater found");
						local_cheating_edges[it->first] ++;
					}
					split_coeff.push_back(((double)it->second)/old_paired_coverage[it->first]);
//		}
			}
			else {
				DEBUG("Zero covered pair info");
			}
		}
		DEBUG("split_edge size " << split_edge.size());
		if (split_edge.size() > 0) {
			pair<VertexId, vector<pair<EdgeId, EdgeId>>> split_pair = new_graph.SplitVertex(v, split_edge, split_coeff);
			res.push_back(split_pair.first);
			if (rc_mode) {
				for( auto it = split_pair.second.begin(); it != split_pair.second.end(); ++it) {
					WrappedSetCoverage(conj_wrap(new_graph, it->second), new_graph.coverage(it->second) * new_graph.length(it->second));
				}
			}
			map<EdgeId, EdgeId> old_to_new_edgeId;
			for(auto it = split_pair.second.begin(); it != split_pair.second.end(); ++it){
				old_to_new_edgeId[it->first] = it->second;
				edge_labels[it->second] = edge_labels[it->first];
				if (cheating_mode) {
					if (local_cheating_edges.find(it->first) != local_cheating_edges.end()) {
						if (local_cheating_edges[it->first] == 0 ) {
							WARN(" 0 copy of edge "<< new_IDs.ReturnIntId(it->first) << " , something wrong");
						} else {
							if (local_cheating_edges[it->first] == 1 ) {
								DEBUG( "cheating OK, no global cheaters needed(but actually added)");
							} else{
								DEBUG( "cheating OK");
							}
							global_cheating_edges.insert(it->second);
						}
					}
				}
			}
			for (size_t j = 0; j < edge_infos.size(); j++){
				if (edge_info_colors[j] == i)
					paired_di_data.ReplaceFirstEdge(edge_infos[j].lp, old_to_new_edgeId[edge_infos[j].lp.first]);
			}
			for(auto it = split_pair.second.begin(); it != split_pair.second.end(); ++it){
				if (new_graph.coverage(it->second) < cfg::get().simp.ec.max_coverage * 0.1) {
				    paired_di_data.DeleteEdgeInfo(it->second);
				    new_graph.DeleteEdge(it->second);
				}
			}
			if (rc_mode)
				BanRCVertex(split_pair.first);
		}
	}
	DEBUG("split finished, deleting vertex");
	new_graph.ForceDeleteVertex(v);
	return res;

}
template<class Graph>
void RepeatResolver<Graph>::BanRCVertex(VertexId v ){
	int id = new_IDs.ReturnIntId(v);
	VertexId rv = conj_wrap(new_graph, v);
	int rc_id = new_IDs.ReturnIntId(rv);
	DEBUG("added vertex " << id << " banning vertex "<< rc_id);
	vector<EdgeId> tmp = new_graph.IncomingEdges(rv);
	for(size_t i = 0; i < tmp.size(); i++)
		global_cheating_edges.insert(tmp[i]);
	TRACE("incoming cheaters added");
	tmp = new_graph.OutgoingEdges(rv);
	for(size_t i = 0; i < tmp.size(); i++)
		global_cheating_edges.insert(tmp[i]);
	TRACE("outgoing cheaters added");
}



template<class Graph>
ConjugateDeBruijnGraph::VertexId RepeatResolver<Graph>::conj_wrap(ConjugateDeBruijnGraph& g, ConjugateDeBruijnGraph::VertexId v){
	return g.conjugate(v);
}

template<class Graph>
NonconjugateDeBruijnGraph::VertexId RepeatResolver<Graph>::conj_wrap(NonconjugateDeBruijnGraph& g, NonconjugateDeBruijnGraph::VertexId v){
	VERIFY(0);
	return v;
}

template<class Graph>
ConjugateDeBruijnGraph::EdgeId RepeatResolver<Graph>::conj_wrap(ConjugateDeBruijnGraph& g, ConjugateDeBruijnGraph::EdgeId e){
	return g.conjugate(e);
}

template<class Graph>
NonconjugateDeBruijnGraph::EdgeId RepeatResolver<Graph>::conj_wrap(NonconjugateDeBruijnGraph& g, NonconjugateDeBruijnGraph::EdgeId e){
	VERIFY(0);
	return e;
}
template<class Graph>
void RepeatResolver<Graph>::WrappedSetCoverage(EdgeId e, int cov){
	if (rc_mode == 0) {
		new_graph.coverage_index().SetCoverage(e, cov);
	} else {
		new_graph.coverage_index().SetCoverage(e, cov);
		EdgeId rc_e = conj_wrap(new_graph, e);
		new_graph.coverage_index().SetCoverage(rc_e, cov);

	}

}


template<class Graph>
map<int, typename Graph::VertexId> RepeatResolver<Graph>::fillVerticesAuto(){
	map<int, typename Graph::VertexId> vertices;
	vertices.clear();
	for (auto v_iter = new_graph.SmartVertexBegin(); !v_iter.IsEnd();
			++v_iter) {
		//			if (vertices.find(new_graph.conjugate(*v_iter)) == vertices.end()) {
		//				if (new_graph.OutgoingEdgeCount(*v_iter) + new_graph.IncomingEdgeCount(*v_iter) > 0)
		vertices.insert(make_pair(10000 - new_IDs.ReturnIntId(*v_iter), *v_iter));
		//			}
	}
	return vertices;
}

template<class Graph>
map<int, typename Graph::VertexId> RepeatResolver<Graph>::fillVerticesComponents(){
	map<int, typename Graph::VertexId> vertices;
	vertices.clear();
	LongEdgesExclusiveSplitter<Graph> splitter(new_graph, cfg::get().ds.IS);

	vector<VertexId> comps;
	if (! splitter.Finished())
		comps = splitter.NextComponent();
	int count = 0;
	while (comps.size() != 0) {
		for(size_t i = 0; i < comps.size(); i++) {
			vertices.insert(make_pair(10000 - new_IDs.ReturnIntId(comps[i]), comps[i]));
			count++;
		}
		if (splitter.Finished())
			break;
		comps = splitter.NextComponent();
	}
	return vertices;

}

template<class Graph>
void RepeatResolver<Graph>::ResolveRepeats(const string& output_folder) {
	//	old_graph = g;
	//	old_index = ind;
	INFO("resolve_repeats started");
	sum_count = 0;
	global_cheating_edges.clear();



	TotalLabelerGraphStruct<Graph> graph_struct_before(old_graph,
			&old_IDs, &old_pos, NULL);
	TotalLabelerGraphStruct<Graph> graph_struct_after(new_graph,
			&new_IDs, &new_pos, NULL);
	TotalLabeler<Graph> TotLabAfter(&graph_struct_after,
			&graph_struct_before);

	for (cheating_mode = 0; cheating_mode < cfg::get().rr.mode; cheating_mode++) {
		INFO(" cheating_mode = " << cheating_mode);
		bool changed = true;
		map<int, VertexId> vertices;

		while (changed) {
			changed = false;
			if (rc_mode)
				vertices = fillVerticesComponents();
			else
				vertices = fillVerticesAuto();
			INFO(
					"Having "<< vertices.size() << " paired vertices, trying to split");
			RealIdGraphLabeler<Graph> IdTrackLabelerAfter(new_graph, new_IDs);
			int GraphCnt = 0;

			omnigraph::WriteSimple(
					new_graph, TotLabAfter,
					output_folder + "resolve_" + ToString(cheating_mode)+"_"+ ToString(GraphCnt) + ".dot",
					"no_repeat_graph");

			for (auto v_iter = vertices.begin(), v_end =
					vertices.end(); v_iter != v_end; ++v_iter) {

				DEBUG(" resolving vertex "<<new_IDs.ReturnIntId(v_iter->second));
				if (rc_mode && deleted_handler.live_vertex.find(v_iter->second) == deleted_handler.live_vertex.end()){
					DEBUG("already deleted");
					continue;
				} else {
					DEBUG("not deleted");
				}
				vector<EdgeId> edgeIds[2];
				int flag = 1;
				edgeIds[0] = new_graph.OutgoingEdges(v_iter->second);
				edgeIds[1] = new_graph.IncomingEdges(v_iter->second);
				for(int i = 0; i < 2; i++) {
					for(size_t j = 0; j < edgeIds[i].size(); j++)
						if (global_cheating_edges.find(edgeIds[i][j]) != global_cheating_edges.end()) {
							flag = 0;
							break;
						}
				}
				if (! flag) {
					DEBUG("Cheaters are near" << new_IDs.ReturnIntId(v_iter->second));
					continue;
				}
				size_t p_size = GenerateVertexPairedInfo(new_graph, paired_di_data, v_iter->second);
				DEBUG("paired info size: " << p_size);

				int tcount;
				if (cheating_mode != 1)
					tcount = RectangleResolveVertex(v_iter->second);
				else
					tcount = CheatingResolveVertex(v_iter->second);
				DEBUG("Vertex "<< v_iter->first<< " resolved to "<< tcount);
				sum_count += tcount;
				if (tcount > 1) {
					GraphCnt++;
					omnigraph::WriteSimple(
						new_graph, TotLabAfter, output_folder + "resolve_" + ToString(cheating_mode)+"_" + ToString(GraphCnt)
								+ ".dot", "no_repeat_graph");
				}
			}
		}
	}INFO("total vert" << sum_count);
	INFO("Converting position labels");

/*	for (auto e_iter = new_graph.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		EdgeId old_edge = edge_labels[*e_iter];
		for (size_t i = 0; i < old_pos.EdgesPositions[old_edge].size(); i++) {
			new_pos.AddEdgePosition(*e_iter,
					old_pos.EdgesPositions[old_edge][i].start(),
					old_pos.EdgesPositions[old_edge][i].end());
		}
	}
*/	//	gvis::WriteSimple(  "repeats_resolved_siiimple.dot", "no_repeat_graph", new_graph);
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

	if (pair_inf.d - new_graph.length(left_id) > 1.3 * cfg::get().ds.IS ) {
		TRACE(
				"PairInfo "<<edge_labels[left_id]<<"("<<new_graph.length(left_id)<<")"<<" "<<right_id<<"("<<old_graph.length(right_id)<<")"<<" "<<pair_inf.d);
//				DEBUG("too far to correct");
		return make_pair(false, pair_inf);
	}

	PairInfo corrected_info = StupidPairInfoCorrectorByOldGraph(new_graph,
			pair_inf);
	TRACE("PairInfo "<<edge_labels[left_id]<<" "<<right_id<<" "<<pair_inf.d<< " corrected into "<<corrected_info.d);
	if (abs(corrected_info.d - pair_inf.d) > MAX_DISTANCE_CORRECTION) {
		TRACE("big correction");
		return make_pair(false, corrected_info);
	}
//	if (corrected_info.d - new_graph.length(left_id) > 130) {
//		DEBUG("too far");
//		return make_pair(false, corrected_info);
//	}
	//todo check correctness. right_id belongs to original graph, not to new_graph.
	if (corrected_info.d + new_graph.length(right_id) < 1/(1.3) * cfg::get().ds.IS) {
		TRACE("too close");
		return make_pair(false, corrected_info);
	}
	TRACE("good");
	return make_pair(true, corrected_info);
//	return make_pair(true, pair_inf);
}

template<class Graph>
size_t RepeatResolver<Graph>::GenerateVertexPairedInfo(Graph &new_graph,
		PairInfoIndexData<EdgeId> &paired_data, VertexId vid) {
	DEBUG("---- Generate vertex paired info for:  " << vid <<" ("<<new_IDs.ReturnIntId(vid) <<") -----------------------------");
	edge_infos.clear();
	local_cheating_edges.clear();
	vector<EdgeId> edgeIds[2];
	edgeIds[0] = new_graph.OutgoingEdges(vid);
	edgeIds[1] = new_graph.IncomingEdges(vid);
	vector<set<EdgeId> > paired_edges;
	DEBUG("out: " << edgeIds[0].size()<< "  in:" << edgeIds[1].size());
	paired_edges.resize(edgeIds[0].size() + edgeIds[1].size());
//TODO:: extract this optional parameter- direction of resolve?
//Or due to r-c structure of original graph it doesn't matter?
	int mult = 1;
	set<EdgeId> right_edges;
	map<VertexId, int> right_vertex_degrees;
	for (int dir = 0; dir < 2; dir++) {
		for (int i = 0, n = edgeIds[dir].size(); i < n; i++) {
			PairInfos tmp = paired_di_data.GetEdgeInfos(edgeIds[dir][i]);
			vector<EdgeInfo> tmp_edge_infos;
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
					if (d * mult >= -1e-8) {

						TRACE("PairInfo: " << new_IDs.ReturnIntId(tmp[j].first)<<" "<<old_IDs.ReturnIntId(edge_labels[tmp[j].first]) << " " << old_IDs.ReturnIntId(tmp[j].second) <<" "<< tmp[j].d);

						pair<bool, PairInfo> correction_result =
								CorrectedAndNotFiltered(new_graph, tmp[j]);
						if (!correction_result.first)
							continue;
//						DEBUG(
//								"PairInfo "<<edge_labels[left_id]<<" "<<right_id<<" "<<d<< " corrected into "<<tmp[j].d<< "weight" << tmp[j].weight);
//						DEBUG(
//								"PairInfo: " << old_IDs.ReturnIntId(edge_labels[tmp[j].first]) << " " << old_IDs.ReturnIntId(tmp[j].second) <<" "<< tmp[j].d);
						EdgeInfo ei(correction_result.second, dir, right_id,
								correction_result.second.d - dif_d);
						int trusted_dist = cfg::get().ds.IS - cfg::get().ds.RL;
						if (cheating_mode == 2 && ((correction_result.second.d - dif_d + old_graph.length(right_id) < trusted_dist - near_vertex) || (correction_result.second.d - dif_d > trusted_dist  + near_vertex))) {
							local_cheating_edges.insert(make_pair(left_id, 0));
							DEBUG("ignored paired_info between " << new_IDs.ReturnIntId(left_id) <<" and " <<old_IDs.ReturnIntId(right_id) <<" with distance " << correction_result.second.d - dif_d);
						} else {
							tmp_edge_infos.push_back(ei);
						//					DEBUG(right_id);
							right_edges.insert(right_id);
						}
					}

				}
			}
			{
				//incompateble with cheating mode
				bool NotOnSelfExist = false;
				for (int j = 0; j < (int)tmp_edge_infos.size(); j++) {
					if ((tmp_edge_infos[j].lp.d != 0)||(edge_labels[tmp_edge_infos[j].lp.first] != tmp_edge_infos[j].lp.second)){
						NotOnSelfExist = true;
						break;
					}
				}
				NotOnSelfExist = false;
				for (int j = 0; j < (int)tmp_edge_infos.size(); j++) {
					edge_infos.push_back(tmp_edge_infos[j]);
				}
			}
		}
	}

	for (int j = 0; j < (int)edge_infos.size(); j++) {
		PairInfo tmp = edge_infos[j].lp;
		DEBUG("Edge infos "<<j<<":"  << new_IDs.ReturnIntId(tmp.first)<<" ("<<old_IDs.ReturnIntId(edge_labels[tmp.first]) << ") -- " << old_IDs.ReturnIntId(tmp.second) <<" "<< tmp.d << " from vertex: "<<edge_infos[j].d<< " weigth "<<tmp.weight);
	}


	return right_edges.size();
}



template<class Graph>
size_t RepeatResolver<Graph>::RectangleResolveVertex(VertexId vid) {
	DEBUG("Rectangle resolve vertex started");
	int size = edge_infos.size();
	//No resolve for cheater-incident vertixes")'
	if (cheating_mode) {
		vector<EdgeId> edgeIds[2];
		edgeIds[0] = new_graph.OutgoingEdges(vid);
		edgeIds[1] = new_graph.IncomingEdges(vid);
		for(int i = 0; i < 2; i++) {
			for(size_t j = 0; j <edgeIds[i].size(); j++ ) {
				if (global_cheating_edges.find(edgeIds[i][j]) != global_cheating_edges.end()) {
					DEBUG ("Can not resolve vertex " << new_IDs.ReturnIntId(vid) <<" because of incident cheater edge "<< new_IDs.ReturnIntId(edgeIds[i][j]));
					return 1;
				}
				if (size == 0) {
					DEBUG ("Can not resolve vertex " << new_IDs.ReturnIntId(vid) <<" because of zero sized info ");
					return 1;
				}
			}
		}

	}
	edge_info_colors.resize(size);
	for (int i = 0; i < size; i++) {
		edge_info_colors[i] = -1;
	}
	vector<vector<int> > neighbours;
	neighbours.resize(size);


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
			if (edge_infos[i].isAdjacent(edge_infos[j], old_graph, new_graph, labels_after)
					&& !edge_infos[j].isAdjacent(edge_infos[i], old_graph,
							new_graph, labels_after))
				WARN(
						"ASSYMETRIC: " << new_IDs.ReturnIntId(edge_infos[i].getEdge()) << " " << new_IDs.ReturnIntId(edge_infos[j].getEdge()));
			if (edge_infos[i].isAdjacent(edge_infos[j], old_graph, new_graph, labels_after)) {
				neighbours[i].push_back(j);
				neighbours[j].push_back(i);
				TRACE(
						old_IDs.ReturnIntId(edge_infos[i].lp.second) <<" " << edge_infos[i].d << " is adjacent "<<old_IDs.ReturnIntId( edge_infos[j].lp.second) <<" " << edge_infos[j].d);

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
	DEBUG("Edge color info " << debruijn::operator<<(oss_, edge_info_colors));
	if (cheating_mode) {
		if (cur_color > 1) {
			DEBUG("cheat_2 resolved vertex " << new_IDs.ReturnIntId(vid));
		} else {
			DEBUG("cheat_2 ignored vertex " << new_IDs.ReturnIntId(vid));
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
				DEBUG("neighbours "<<first_ind<<" + "<<counts[1]<<"  "<<second_ind);
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
	DEBUG("Colours " << debruijn::operator<<(oss_, cheater_colors));

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
