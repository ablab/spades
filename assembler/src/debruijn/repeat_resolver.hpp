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
#include <algorithm>
#include "logging.hpp"
#include "paired_info.hpp"
#include "config.hpp"
#include "omni_utils.hpp"

LOGGER("d.repeat_resolver");
namespace de_bruijn {
using omnigraph::SmartVertexIterator;

template<class Graph>
class RepeatResolver {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	typedef de_bruijn::SmartVertexIterator<Graph> VertexIter;
	typedef omnigraph::SmartEdgeIterator<Graph> EdgeIter;
	typedef de_bruijn::PairedInfoIndex<Graph> PairedInfoIndex;
	typedef typename PairedInfoIndex::PairInfo PairedInfo;
	typedef vector<PairedInfo> PairInfos;

	typedef map<VertexId,set<EdgeId> > NewVertexMap;
	typedef map <VertexId, set<VertexId> > VertexIdMap;
public:
	RepeatResolver(int leap = 0) : leap_(leap), new_graph(K), new_index(new_graph){


		assert(leap >= 0 && leap < 100);
	}
	Graph ResolveRepeats(Graph &g, PairedInfoIndex &ind);

private:
	int leap_;
	size_t ResolveVertex(Graph &g, PairedInfoIndex &ind, VertexId vid);
	void ResolveEdge(EdgeId eid);
	void dfs(vector<vector<int> > &edge_list, vector<int> &colors, int cur_vert, int cur_color);
	VertexIdMap vid_map;
	NewVertexMap new_map;
	Graph new_graph;
	int sum_count;
//	PairedInfoIndex old_index;
	PairedInfoIndex new_index;

//	Graph old_graph;
};

template<class Graph>
Graph RepeatResolver<Graph>::ResolveRepeats(Graph &g, PairedInfoIndex &ind){
//	old_graph = g;
//	old_index = ind;
	INFO("resolve_repeats started");
	int count = 0;
	sum_count = 0;
	bool  changed = true;
	set<VertexId> vertices;
	while (changed) {
		changed = false;
		vertices.clear();
		for(VertexIter v_iter = g.SmartVertexBegin(), end = g.SmartVertexEnd(); v_iter != end; ++v_iter) {
			if (vertices.find(g.conjugate(*v_iter)) == vertices.end())
			{
				vertices.insert(*v_iter);
			}
		}
		INFO("Having "<< vertices.size() << "paired vertices, trying to split");
		for(auto v_iter = vertices.begin(), v_end = vertices.end(); v_iter != v_end; ++v_iter) {
//		vector<typename PairedInfoIndex::PairInfo> tmp = old_index.getEdgeInfos(*e_iter);
//		INFO("Parsing vertex "<< old_graph.VertexNucls(*v_iter));
			if (ResolveVertex(g, ind, *v_iter)  > 1) {
// 				changed = true;
			}
			TRACE("Vertex "<< count << " resolved");
		}
	}

	//	for(EdgeIter e_iter = old_graph.SmartEdgeBegin(), end = old_graph.SmartEdgeEnd(); e_iter != end; ++e_iter) {
//		PairInfos tmp = old_index.GetEdgeInfo(*e_iter);
		//		ResolveEdge(*e_iter);G
//	}
	INFO("total vert" << sum_count);

	return new_graph;
}

template<class Graph>
void RepeatResolver<Graph>::dfs(vector<vector<int> > &edge_list, vector<int> &colors, int cur_vert, int cur_color){
	colors[cur_vert] = cur_color;
	DEBUG("dfs-ing, vert num" << cur_vert << " " << cur_color);
	for(int i = 0, sz = edge_list[cur_vert].size(); i <sz; i ++) {
		if (colors[edge_list[cur_vert][i]] > 0) {
			if (colors[edge_list[cur_vert][i]] != cur_color) {
				ERROR("error in dfs, neighbour to " << edge_list[cur_vert][i]);
			}
		} else {
			dfs(edge_list, colors, edge_list[cur_vert][i], cur_color);
		}
	}
}
template<class Graph>
size_t RepeatResolver<Graph>::ResolveVertex( Graph &g, PairedInfoIndex &ind, VertexId vid){
	DEBUG("Parsing vertex " << vid);
	DEBUG(g.conjugate(vid));
	vector<EdgeId> edgeIds[2];
	edgeIds[0] = g.OutgoingEdges(vid);
	edgeIds[1] = g.IncomingEdges(vid);
	vector<set<EdgeId> > paired_edges;
	DEBUG(edgeIds[0].size()<< "  " << edgeIds[1].size());
	paired_edges.resize(edgeIds[0].size() + edgeIds[1].size());

	map<EdgeId, set<EdgeId> > right_to_left;
	map<EdgeId, int> right_set;
	vector<EdgeId> right_vector[2];
	map<EdgeId, set<int> > left_colors[2];
	vector<int> colors[2];
	for (int mult = -1; mult < 3; mult += 2) {
		int cur_ind = (mult + 1) / 2;
		int cur_id = 0;

		right_set.clear();
		right_vector[cur_ind].clear();
		right_to_left.clear();
		for (int dir = 0; dir < 2; dir++) {
			for (int i = 0, n = edgeIds[dir].size(); i < n; i ++) {
				PairInfos tmp = ind.GetEdgeInfo(edgeIds[dir][i]);
				for (int j = 0, sz = tmp.size(); j < sz; j++) {
					EdgeId right_id = tmp[j].second();
					EdgeId left_id = tmp[j].first();
					if (tmp[j].d() * mult > 0) {
						if (right_to_left.find(right_id) != right_to_left.end())
							right_to_left[right_id].insert(left_id);
						else {
							set<EdgeId> tmp_set;
							tmp_set.insert(left_id);
							right_to_left.insert(make_pair(right_id, tmp_set));
							right_set.insert(make_pair(right_id, cur_id));
							right_vector[cur_ind].push_back(right_id);
							cur_id ++;
						}

					}
				}
		//		old_index.getEdgeInfos(inEdgeIds[i]);
			}
		}



		size_t right_edge_count = right_set.size();
		vector<vector<int> > edge_list(right_edge_count);
		DEBUG("Total: " << right_edge_count << "edges");
		LOG_ASSERT(right_edge_count == right_vector[cur_ind].size(), "Size mismatch");
		colors[cur_ind].resize(right_edge_count);
		for(size_t i = 0; i < right_edge_count; i++)
			colors[cur_ind][i] = 0;
		for(size_t i = 0; i < right_edge_count; i++) {
	//TODO Add option to "jump" - use not only direct neighbours(parameter leap in constructor)
			vector<EdgeId> neighbours = g.NeighbouringEdges(right_vector[cur_ind][i]);
			DEBUG("neighbours to "<<i<<"(" <<right_vector[cur_ind][i]<<  "): " << neighbours.size())


			for(int j = 0, sz = neighbours.size(); j < sz; j++){
				if (right_set.find(neighbours[j]) != right_set.end()) {
					edge_list[i].push_back(right_set[neighbours[j]]);
					DEBUG (right_set[neighbours[j]]<<"  "<<neighbours[j]);
				}
			}
		}
		int cur_color = 0;
		DEBUG("dfs started");
		for(size_t i = 0; i < right_edge_count; i++) {
			if (colors[cur_ind][i] == 0) {
				cur_color++;
				dfs(edge_list, colors[cur_ind], i, cur_color);
			}
		}
		for(size_t i = 0; i < right_edge_count; i++) {
			for(typename set<EdgeId>::iterator e_it = right_to_left[right_vector[cur_ind][i]].begin() , end_it = right_to_left[right_vector[cur_ind][i]].end(); e_it != end_it; e_it ++){
				if (left_colors[cur_ind].find(*e_it) == left_colors[cur_ind].end()) {
					set<int> tmp;
					left_colors[cur_ind].insert(make_pair(*e_it, tmp));
				}
				left_colors[cur_ind][*e_it].insert(colors[cur_ind][i]);
			}
		}
	}
	map<pair<int, int> ,int> color_pairs;
	int new_colors_count = 0;
	for(auto iter = left_colors[0].begin(), end_iter = left_colors[0].end();iter != end_iter; iter++ ) {
		typename map<EdgeId, set<int> >::iterator sec_iter;
		sec_iter = left_colors[1].find(iter->first);
		if (sec_iter == left_colors[1].end()) {
//TODO:: accuratly csheck this boundary cases
			WARN("an edge " << iter->first << " has paired info in only one direction");
			continue;
		}

		for(auto first_c_iter = iter->second.begin(),first_c_end_iter = iter->second.end(); first_c_iter != first_c_end_iter; first_c_iter ++ ){
			for(auto second_c_iter = sec_iter->second.begin(),second_c_end_iter = sec_iter->second.end(); second_c_iter != second_c_end_iter; second_c_iter ++ ){
				pair<int, int> tmp = make_pair(*first_c_iter, *second_c_iter);
				if (color_pairs.find(tmp) == color_pairs.end()) {
					color_pairs.insert(make_pair(tmp, new_colors_count));
					new_colors_count ++;
				}
			}
		}
	}
	if (color_pairs.size() > 1) {
		TRACE("vertex "<<vid<< " splitted to " << color_pairs.size());
	}

	if (new_colors_count > 1) {

		vector<VertexId> new_verts;
		new_verts.resize(new_colors_count);
		for(int i = 0; i < new_colors_count ; i++){
			new_verts[i] = g.AddVertex();
		}
		for (auto e_iter0 = left_colors[0].begin(), end_iter = left_colors[0].end(); e_iter0 != end_iter; e_iter0++) {
			typename map<EdgeId, set<int> >::iterator e_iter1 = left_colors[1].find(e_iter0->first);
//TODO: no left coverage must be resolved

			if (e_iter1 != left_colors[1].end() ) {
				for (auto color0_iter = e_iter0->second.begin(); color0_iter != e_iter0->second.end(); color0_iter++){
					for (auto color1_iter = e_iter1->second.begin(); color1_iter != e_iter1->second.end(); color1_iter++){
						VertexId new_f[2];
						new_f[0] = g.EdgeStart(e_iter0->first);
						new_f[1] = g.EdgeEnd(e_iter0->first);
						for(int ii = 0; ii < 2; ii++) {
							if (new_f[ii] == vid)
								new_f[ii] = new_verts[color_pairs[make_pair(*color0_iter,*color1_iter)]];
						}
						g.AddEdge(new_f[0], new_f[1], g.EdgeNucls(e_iter0->first));
					}

				}
			}
		}
		for(int ii = 0; ii < 2; ii++)
			for(auto e_iter = edgeIds[ii].begin(), end_iter = edgeIds[ii].end(); e_iter != end_iter; e_iter++) {
				g.DeleteEdge(*e_iter);
			}
		g.DeleteVertex(vid);
	}
	sum_count += color_pairs.size();
	return color_pairs.size();

}

}
#endif /* REPEAT_RESOLVER_HPP_ */
