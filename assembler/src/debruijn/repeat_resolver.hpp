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

namespace debruijn_graph {

using omnigraph::SmartVertexIterator;
using omnigraph::Compressor;


template<class Graph>
class RepeatResolver {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	typedef SmartVertexIterator<Graph> VertexIter;
	typedef omnigraph::SmartEdgeIterator<Graph> EdgeIter;
	typedef PairedInfoIndex<Graph> PIIndex;
	typedef vector<PairInfo> PairInfos;

	typedef map<VertexId,set<EdgeId> > NewVertexMap;
	typedef map <VertexId, set<VertexId> > VertexIdMap;

public:

	class EdgeInfo {

	public:
		static const int MAXD = 100;


		EdgeInfo(const PairInfo &lp_, const int dir_ , const EdgeId edge_, const int d_) : lp(lp_), dir(dir_), edge(edge_), d(d_){

		}
		inline EdgeId getEdge(){
			return edge;
		}

		bool isClose(int a, int b){
			return (abs(a - b) < MAXD);
		}
		bool isAdjacent(EdgeInfo other_info, Graph &old_graph) {
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
			if ((v_e == other_v_s && isClose(d + len, other_d))||( v_s == other_v_e && isClose(d, other_d + other_len)) ||( other_edge == edge  && isClose(d, other_d))) {
//				DEBUG("ADJACENT!");
				return true;
			}
			else {
//			DEBUG("not adjacent");
				return false;
			}
		}

		inline int getDistance(){
			return d;
		}
	public:
		PairInfo lp;
		int dir;
	private:
		EdgeId edge;
		int d;

	};




	RepeatResolver(Graph &old_graph_, int leap, PIIndex &ind, Graph &new_graph_ ) : leap_(leap), new_graph(new_graph_), old_graph(old_graph_){
		unordered_map<VertexId, VertexId> old_to_new;
		unordered_map<EdgeId, EdgeId> old_to_new_edge;
		size_t paired_size = 0;
		set<VertexId> vertices;
		vertices.clear();
		real_vertices.clear();
		set<EdgeId> edges;
		edges.clear();
		for(VertexIter v_iter = old_graph.SmartVertexBegin(); !v_iter.IsEnd(); ++v_iter) {
	//		if (vertices.find(old_graph.conjugate(*v_iter)) == vertices.end())
			{
				vertices.insert(*v_iter);
			}
		}
		for(EdgeIter e_iter = old_graph.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		//	if (edges.find(old_graph.conjugate(*e_iter)) == edges.end())
			{
				edges.insert(*e_iter);
			}
		}
		for(auto v_iter = vertices.begin();v_iter != vertices.end(); ++v_iter) {
			size_t degree = old_graph.IncomingEdgeCount(*v_iter) + old_graph.OutgoingEdgeCount(*v_iter);
		    if (degree > 0)
			{
				VertexId new_vertex = new_graph.AddVertex();
				real_vertices.insert(new_vertex);
				DEBUG("Added vertex" << new_vertex <<" " << new_graph.conjugate(new_vertex));
				vertex_labels[new_vertex] = *v_iter;
				old_to_new[*v_iter] = new_vertex;
			//	old_to_new[old_graph.conjugate(*v_iter)] = new_graph.conjugate(new_vertex);
			}

		}
		for(auto e_iter = edges.begin(); e_iter != edges.end(); ++e_iter) {
			DEBUG("Adding edge from " << old_to_new[old_graph.EdgeStart(*e_iter)] <<" to " << old_to_new[old_graph.EdgeEnd(*e_iter)]);
//			DEBUG("Setting coverage to edge length " << old_graph.length(*e_iter) << "  cov: " << old_graph.coverage(*e_iter));
			EdgeId new_edge = new_graph.AddEdge(old_to_new[old_graph.EdgeStart(*e_iter)], old_to_new[old_graph.EdgeEnd(*e_iter)], old_graph.EdgeNucls(*e_iter), old_graph.coverage(*e_iter) * old_graph.length(*e_iter));
			edge_labels[new_edge] = *e_iter;
			DEBUG("Adding edge " << new_edge<< " from" << *e_iter);
			old_to_new_edge[*e_iter] = new_edge;
//			old_to_new_edge[old_graph.conjugate(*e_iter)] = new_graph.conjugate(new_edge);
//			PairInfos tmp = ind.GetEdgeInfo(edgeIds[dir][i]);

		}
		old_to_new.clear();
		for (auto p_iter = ind.begin(), p_end_iter = ind.end(); p_iter != p_end_iter; ++p_iter){
			PairInfos pi = *p_iter;
			paired_size += pi.size();
			for (size_t j = 0; j < pi.size(); j++) {
				if (old_to_new_edge.find(pi[j].first) != old_to_new_edge.end() && old_to_new_edge.find(pi[j].second) != old_to_new_edge.end()) {
					TRACE("Adding pair " << pi[j].first<<"  " <<old_to_new_edge[pi[j].first] << "  " <<pi[j].second);
					PairInfo *tmp = new PairInfo(old_to_new_edge[pi[j].first], pi[j].second, pi[j].d, pi[j].weight);
					paired_di_data.AddPairInfo(*tmp, 0);
				} else {
					WARN("Paired Info with deleted edge! " << pi[j].first<<"  "  <<pi[j].second);
				}
			}
		}
		INFO("paired info size: "<<paired_size);
		assert(leap >= 0 && leap < 100);
	}
	void ResolveRepeats();

private:
	int leap_;
	size_t RectangleResolveVertex(VertexId vid);

	size_t GenerateVertexPairedInfo(Graph &g, PairInfoIndexData &ind, VertexId vid);
	vector<typename Graph::VertexId> MultiSplit(VertexId v);

	void ResolveEdge(EdgeId eid);
	void dfs(vector<vector<int> > &edge_list, vector<int> &colors, int cur_vert, int cur_color);
	VertexIdMap vid_map;
	NewVertexMap new_map;
	Graph &new_graph;
	Graph &old_graph;
	vector<int> edge_info_colors;
	vector<EdgeInfo> edge_infos;
	PairInfoIndexData paired_di_data;
	unordered_map<VertexId, VertexId> vertex_labels;
	unordered_map<EdgeId, EdgeId> edge_labels;
	set<VertexId> real_vertices;

	int sum_count;

private:
	DECL_LOGGER("RepeatResolver")
};
template<class Graph>
vector<typename Graph::VertexId> RepeatResolver<Graph>::MultiSplit(VertexId v){
	int k = 0;
	for(size_t i = 0; i < edge_info_colors.size(); i++)
		if (edge_info_colors[i] >= k)
			k ++;
	DEBUG("splitting to "<< k <<" parts");
	vector<VertexId> res;
	res.resize(k);
	vector<EdgeId> edgeIds[2];
//TODO: fix labels
	edgeIds[0] = new_graph.OutgoingEdges(v);
	edgeIds[1] = new_graph.IncomingEdges(v);
	vector<unordered_map<EdgeId, EdgeId> > new_edges(k);
	unordered_map<EdgeId, int> old_paired_coverage;
	unordered_map<EdgeId, int> new_paired_coverage;
	//Remember-because of loops there can be same edges in edgeIds[0] and edgeIds[1]
	for(int i = 0; i < k ; i++) {
		res[i] = new_graph.AddVertex();
	}
	for(size_t i = 0; i < edge_info_colors.size(); i++) {
		EdgeId le = edge_infos[i].lp.first;
		if (old_paired_coverage.find(le) == old_paired_coverage.end())
			old_paired_coverage[le] = 0;
		old_paired_coverage[le]++;
		TRACE("replacing edge " << le<<" with label " << edge_labels[le] << " "<< (new_edges[edge_info_colors[i]].find(le) == new_edges[edge_info_colors[i]].end()));

		EdgeId res_edge = NULL;
		if (edge_infos[i].dir == 0 ) {
			if (new_graph.EdgeStart(le) != v) {
				ERROR("Non incident edge!!!" << new_graph.EdgeStart(le) <<" instead of "<< v);
			} else {

				if (new_edges[edge_info_colors[i]].find(le) == new_edges[edge_info_colors[i]].end())
					res_edge = new_graph.AddEdge(res[edge_info_colors[i]], new_graph.EdgeEnd(le), new_graph.EdgeNucls(le), 0);
			}
		} else {
			if (new_graph.EdgeEnd(le) != v) {
				ERROR("Non incident edge!!!" << new_graph.EdgeEnd(le) <<" instead of "<< v);
			} else {
				if (new_edges[edge_info_colors[i]].find(le) == new_edges[edge_info_colors[i]].end())
					res_edge = new_graph.AddEdge( new_graph.EdgeStart(le),res[edge_info_colors[i]], new_graph.EdgeNucls(le), 0);
			}
		}
		TRACE("replaced");
		if (res_edge != NULL) {
			new_edges[edge_info_colors[i]].insert(make_pair(le, res_edge));
			edge_labels[res_edge] = edge_labels[le];
			TRACE("before replace first Edge");
			paired_di_data.ReplaceFirstEdge(edge_infos[i].lp, res_edge);
			new_paired_coverage[new_edges[edge_info_colors[i]][le]] = 0;
		}
		else {
			paired_di_data.ReplaceFirstEdge(edge_infos[i].lp, new_edges[edge_info_colors[i]][le]);
		}
		new_paired_coverage[new_edges[edge_info_colors[i]][le]] ++;
	}
	for(int i = 0; i < k; i ++) {
		for (auto edge_iter = new_edges[i].begin(); edge_iter != new_edges[i].end(); edge_iter ++) {
			DEBUG("setting coverage, from edgeid "<< edge_iter->first<<" length "<<new_graph.length(edge_iter->first) <<"   "<< new_graph.coverage(edge_iter->first) <<" taking "<< new_paired_coverage[edge_iter->second] <<"/"<< old_paired_coverage[edge_iter->first]);
			new_graph.SetCoverage(edge_iter->second, new_graph.length(edge_iter->first) * new_graph.coverage(edge_iter->first) * new_paired_coverage[edge_iter->second]/ old_paired_coverage[edge_iter->first]);
		}
	}

	new_graph.ForceDeleteVertex(v);
	return res;

}

template<class Graph>
void RepeatResolver<Graph>::ResolveRepeats(){
//	old_graph = g;
//	old_index = ind;
	INFO("resolve_repeats started");
	sum_count = 0;
	bool  changed = true;
	set<VertexId> vertices;
	while (changed) {
		changed = false;
		vertices.clear();
		for(VertexIter v_iter = new_graph.SmartVertexBegin(); !v_iter.IsEnd(); ++v_iter) {
			if (vertices.find(new_graph.conjugate(*v_iter)) == vertices.end())
			{
				vertices.insert(*v_iter);
			}
		}
		INFO("Having "<< vertices.size() << "paired vertices, trying to split");
		for(auto v_iter = real_vertices.begin(), v_end = real_vertices.end(); v_iter != v_end; ++v_iter) {
			DEBUG(" resolving vertex"<<*v_iter<<" "<<GenerateVertexPairedInfo(new_graph, paired_di_data, *v_iter));
			int tcount = RectangleResolveVertex(*v_iter);
			DEBUG("Vertex "<< *v_iter<< " resolved to "<< tcount);
			sum_count += tcount;
		}
	}
	INFO("total vert" << sum_count);
 //	gvis::WriteSimple(  "repeats_resolved_siiimple.dot", "no_repeat_graph", new_graph);
}

template<class Graph>
void RepeatResolver<Graph>::dfs(vector<vector<int> > &edge_list, vector<int> &colors, int cur_vert, int cur_color){
	colors[cur_vert] = cur_color;
//	DEBUG("dfs-ing, vert num" << cur_vert << " " << cur_color);
	for(int i = 0, sz = edge_list[cur_vert].size(); i <sz; i ++) {
		if (colors[edge_list[cur_vert][i]] >  -1) {
			if (colors[edge_list[cur_vert][i]] != cur_color) {
				ERROR("error in dfs, neighbour to " << edge_list[cur_vert][i]);
			}
		} else {
			if (i != cur_vert) {
				dfs(edge_list, colors, edge_list[cur_vert][i], cur_color);
			}
		}
	}
}

template<class Graph>
size_t RepeatResolver<Graph>::GenerateVertexPairedInfo( Graph &new_graph, PairInfoIndexData &paired_data, VertexId vid){
	DEBUG("Generate vertex paired info for:  " << vid);
//	DEBUG(new_graph.conjugate(vid));
	edge_infos.clear();
	vector<EdgeId> edgeIds[2];
	edgeIds[0] = new_graph.OutgoingEdges(vid);
	edgeIds[1] = new_graph.IncomingEdges(vid);
	vector<set<EdgeId> > paired_edges;
	DEBUG(edgeIds[0].size()<< "  " << edgeIds[1].size());
	paired_edges.resize(edgeIds[0].size() + edgeIds[1].size());

	int mult = 1;
	set<EdgeId> neighbours;
	for (int dir = 0; dir < 2; dir++) {
		for (int i = 0, n = edgeIds[dir].size(); i < n; i ++) {
			PairInfos tmp = paired_di_data.GetEdgeInfos(edgeIds[dir][i]);
			for (int j = 0, sz = tmp.size(); j < sz; j++) {
				EdgeId right_id = tmp[j].second;
				EdgeId left_id = tmp[j].first;
				int d = tmp[j].d;
//				int w = tmp[j].weight;
//				if (w < 10) continue;
				int new_d = d;
				if (dir == 1)
					new_d -= new_graph.length(left_id);
				if (d * mult > 0) {
					EdgeInfo ei(tmp[j], dir, right_id, new_d);
					edge_infos.push_back(ei);
//					DEBUG(right_id);
					neighbours.insert(right_id);
				}
			}
		}
	}
	return neighbours.size();
}

template<class Graph>
size_t RepeatResolver<Graph>::RectangleResolveVertex(VertexId vid){
	DEBUG("Rectangle resolve vertex started");
	int size = edge_infos.size();
	edge_info_colors.resize(size);
	for(int i = 0; i < size; i++) {
		edge_info_colors[i] = -1;
	}
	vector<vector<int> > neighbours;
	neighbours.resize(size);
	DEBUG("constructing edge-set");
	set<EdgeId> edges;
	edges.clear();

	for(EdgeIter e_iter = old_graph.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
		{
			edges.insert(*e_iter);
		}
	}
	DEBUG("checking edge-infos ");
	for(int i = 0; i < size; i++) {
		if(edges.find(edge_infos[i].getEdge()) == edges.end()) {
			ERROR("fake edge: "<< edge_infos[i].getEdge());
		}
	}
	for(int i = 0; i < size; i++) {
//		DEBUG("Filling " << i << " element");
		if(edges.find(edge_infos[i].getEdge()) == edges.end()) {
			ERROR("fake edge");
		}
		neighbours[i].resize(0);
		for(int j = 0; j < size; j++){

			if (edge_infos[i].isAdjacent(edge_infos[j], old_graph)){
				neighbours[i].push_back(j);
			}
		}
	}
	int cur_color = 0;
	DEBUG("dfs started");
	for(int i = 0; i < size; i ++){
		if (edge_info_colors[i] == -1) {
			dfs(neighbours, edge_info_colors, i, cur_color);
			cur_color++;
		}
	}
	MultiSplit(vid);
	return cur_color;
}


}

#endif /* REPEAT_RESOLVER_HPP_ */
