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
LOGGER("d.repeat_resolver");


namespace de_bruijn {
template<class Graph>
class RepeatResolver {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	typedef de_bruijn::SmartVertexIterator<Graph> VertexIter;
	typedef de_bruijn::SmartEdgeIterator<Graph> EdgeIter;
	typedef de_bruijn::PairedInfoIndex<Graph> PairedInfoIndex;
	typedef typename PairedInfoIndex::PairInfo PairedInfo;
	typedef vector<PairedInfo> PairInfos;

	typedef map<VertexId,set<EdgeId> > NewVertexMap;
	typedef map <VertexId, set<VertexId> > VertexIdMap;
public:
	RepeatResolver(int leap = 0, Graph g = NULL, PairedInfoIndex old_index_ = NULL) : leap_(leap), new_graph(K), old_graph(g), new_index(new_graph), old_index(old_index_){


		assert(leap >= 0 && leap < 100);
	}
	Graph ResolveRepeats(Graph &g, PairedInfoIndex &ind);

private:
	int leap_;
	void ResolveVertex(Graph &g, PairedInfoIndex &ind, VertexId vid);
	void ResolveEdge(EdgeId eid);
	VertexIdMap vid_map;
	NewVertexMap new_map;
	Graph new_graph;
	PairedInfoIndex old_index;
	PairedInfoIndex new_index;

	Graph old_graph;
};

template<class Graph>
Graph RepeatResolver<Graph>::ResolveRepeats(Graph &g, PairedInfoIndex &ind){
//	old_graph = g;
//	old_index = ind;
	INFO("resolve_repeats started");
	for(VertexIter v_iter = g.SmartVertexBegin(), end = g.SmartVertexEnd(); v_iter != end; ++v_iter) {
//		vector<typename PairedInfoIndex::PairInfo> tmp = old_index.getEdgeInfos(*e_iter);
		INFO("smartiterators ");
//		INFO("Parsing vertex "<< old_graph.VertexNucls(*v_iter));

		ResolveVertex(g, ind, *v_iter);
	}
//	for(EdgeIter e_iter = old_graph.SmartEdgeBegin(), end = old_graph.SmartEdgeEnd(); e_iter != end; ++e_iter) {
//		PairInfos tmp = old_index.GetEdgeInfo(*e_iter);
		//		ResolveEdge(*e_iter);G
//	}
	return new_graph;
}

template<class Graph>
void RepeatResolver<Graph>::ResolveVertex( Graph &g, PairedInfoIndex &ind, VertexId vid){
	INFO("Parsing vertex ");
	INFO("with seq " <<g.VertexNucls(vid));
	vector<EdgeId> edgeIds[2];
	DEBUG("before outgoing edges");
	edgeIds[0] = g.OutgoingEdges(vid);
	edgeIds[1] = g.IncomingEdges(vid);
	vector<set<EdgeId> > paired_edges;
	paired_edges.resize(edgeIds[0].size() + edgeIds[1].size());
	DEBUG("std inited");
	int i;
	unsigned int j;
	i = j;
	map<EdgeId, set<EdgeId>> right_to_left;
	map<EdgeId, int> right_set;
	vector<EdgeId> right_vector;
	int cur_id = 0;
	for (int dir = 0; dir < 2; dir++) {
		for (int i = 0, n = edgeIds[dir].size(); i < n; i ++) {
			DEBUG("edge " << dir <<" "<<i);
			PairInfos tmp = ind.GetEdgeInfo(edgeIds[dir][i]);
			for (int j = 0, sz = tmp.size(); j < sz; j++) {
				EdgeId right_id = tmp[j].second();
				EdgeId left_id = tmp[j].first();
				if (right_to_left.find(right_id) != right_to_left.end())
					right_to_left[right_id].insert(left_id);
				else {
					set<EdgeId> tmp_set;
					tmp_set.insert(left_id);
					right_to_left.insert(make_pair(right_id, tmp_set));
					right_set.insert(make_pair(right_id, cur_id));
					right_vector.push_back(right_id);
					cur_id ++;
				}
			}
	//		old_index.getEdgeInfos(inEdgeIds[i]);
		}
	}
	int right_edge_count = right_set.size();
	vector<vector<EdgeId> > edge_list(right_edge_count);
	LOG_ASSERT(right_edge_count == right_vector.size(), "Size mismatch");
	vector<int> colors(right_edge_count);
	for(int i = 0; i < right_edge_count; i++)
		colors[i] = 0;
	for(int i = 0; i < right_edge_count; i++) {
//TODO Add option to "jump" - use not only direct neighbours(parameter leap in constructor)
		vector<EdgeId> neighbours = g.NeighbouringEdges(right_vector[i]);
		for(int j = 0, sz = neighbours.size(); j < sz; j++){
			if (right_set.find(neighbours[j]) != right_set.end()) {
				edge_list[i].push_back(neighbours[j]);
			}
		}
	}
}

}
#endif /* REPEAT_RESOLVER_HPP_ */


//void dfs (int **table, int color, int * leftcolor, int* rightcolor,  int pos) {
//	leftcolor[pos] = color;
//	forn(j, MAX_DEGREE) {
//		if (table[pos][j] && !(rightcolor[j])) {
//			rightcolor[j] = color;
//			forn(i, MAX_DEGREE)
//				if (table[i][j] && i != pos && !leftcolor[i]) {
//					leftcolor[i] = color;
//					dfs(table,color,leftcolor,rightcolor,i);
//				}
//		}
//	}
//}
//void doSplit(PairedGraph &graph, edgePairsMap &EdgePairs) {
//	int *table[MAX_DEGREE];
//	forn(i, MAX_DEGREE)
//		table[i] = new int[MAX_DEGREE];
//	int leftcolor [MAX_DEGREE];
//	int rightcolor [MAX_DEGREE];
//	map<int, int> edgeIds;
//	map<int, int> idEdges;
//	for(edgePairsMap::iterator iter = EdgePairs.begin(); iter != EdgePairs.end(); iter++) {
//		int curVId = iter->first;
//		INFO(" splitting vertex" << curVId);
//		int len = iter->second.size();
//		forn(i, MAX_DEGREE) {
//			forn(j, MAX_DEGREE)
//				table[i][j] = 0;
//			leftcolor[i] = 0;
//			rightcolor[i] = 0;
//		}
//		cerr<<"Vertex "<<curVId<<" in degree "<<graph.degree(curVId,LEFT)<<" out degree "<<graph.degree(curVId,RIGHT)<<" edge pairs:"<<endl;
//		forn(i, len){
//			cerr<<iter->second[i].first<<" "<<iter->second[i].second<<endl;
//		}
//
//		edgeIds.clear();
//		idEdges.clear();
//		int leftId = 1;
//		int rightId = 1;
//		forn(i, len) {
//			if (edgeIds.find(iter->second[i].first) == edgeIds.end()){
//				edgeIds[iter->second[i].first] = leftId;
//				idEdges[leftId] = iter->second[i].first;
//				leftId ++;
//			}
//			if (edgeIds.find(-iter->second[i].second) == edgeIds.end()){
//				edgeIds[-iter->second[i].second] = rightId;
//				idEdges[-rightId] = iter->second[i].second;
//				rightId ++;
//			}
//		//	leftGlobalIds[leftId] = mp(iter->second[i].first,
//			table[edgeIds[iter->second[i].first]][edgeIds[-iter->second[i].second]] = 1;
//		}
//		int in_degree = graph.degree(curVId, LEFT);
//		int out_degree = graph.degree(curVId, RIGHT);
//		forn(i, in_degree) {
//			if (edgeIds.find(graph.leftEdge(curVId, i)->EdgeId) == edgeIds.end()){
//				edgeIds[graph.leftEdge(curVId, i)->EdgeId] = leftId;
//				idEdges[leftId] = graph.leftEdge(curVId, i)->EdgeId;
//				leftId++;
//
////				EdgePairs.push_back(mp(graph.leftEdge(curVId, i)->EdgeId,0));
//			}
//		}
//		forn(i, out_degree) {
//			if (edgeIds.find(-graph.rightEdge(curVId, i)->EdgeId) == edgeIds.end()){
//				edgeIds[-graph.rightEdge(curVId, i)->EdgeId] = rightId;
//				idEdges[-rightId] = graph.rightEdge(curVId, i)->EdgeId;
//				rightId++;
//			}
//		}
//		int color = 0;
//		for(int i = 1; i < leftId; i++) {
//			if(!leftcolor[i]) {
//
//				color++;
//				dfs(table, color, leftcolor, rightcolor, i);
//			}
//		}
//		for(int i =1; i < rightId; i++) {
//			if(!rightcolor[i]) {
//				color++;
//				rightcolor[i] = color;
//			}
//		}
//		cerr<<"left colors: "<<endl;
//		forn(i, leftId) {
//			cerr<<leftcolor[i]<<" ";
//		}
//		cerr<<"\n right colors: "<<endl;
//		forn(i, rightId) {
//			cerr<<rightcolor[i]<<" ";
//		}
//		cerr<<"table: "<<endl;
//		forn(i, leftId) {
//			forn(j, rightId) {
//			cerr<<table[i][j]<<" ";
//			}
//
//			cerr<<endl;
//		}
//
//
//
//		int tmpVertCount = graph.VertexCount;
//		cerr << "idEdges" << endl;
//		for(map<int, int>::iterator it = idEdges.begin(); it != idEdges.end(); it++)
//			cerr<< it->first << " " << it->second << endl;
//		for(int cur_color = 2;cur_color <= color; cur_color++) {
//
//			int tmp = graph.addVertex(graph.VertexCount + cur_color - 2);
//			INFO("added vertex" << tmp);
//		}
//		forn(i, leftId) {
//			if (i && leftcolor[i] >= 2)
//				graph.addEdgeVertexAdjacency(tmpVertCount + leftcolor[i] - 2, graph.longEdges[idEdges[i]], LEFT);
//		}
//		forn(i, rightId) {
//			if (i && rightcolor[i] >= 2)
//				graph.addEdgeVertexAdjacency(tmpVertCount + rightcolor[i] - 2, graph.longEdges[idEdges[-i]], RIGHT);
//		}
//		//graph.VertexCount += color - 2;
//	}
//	graph.recreateVerticesInfo(graph.VertexCount, graph.longEdges);
//
//}
//
//void SplitByLowers(PairedGraph &graph){
//	edgePairsMap EdgePairs;
//	forn(CurVert, graph.VertexCount){
//		vector<pair<int,int>> tmpVect;
//		forn (edgeI, graph.degree(CurVert, LEFT)){
//			Edge* leftEdge = graph.neighbourEdge(CurVert,edgeI,LEFT);
//			forn (edgeJ, graph.degree(CurVert, RIGHT)){
//				Edge* rightEdge = graph.neighbourEdge(CurVert,edgeJ,RIGHT);
//	//			cerr<<"Check: "<<CurVert<<" Edge pair "<<leftEdge->EdgeId<<" "<<rightEdge->EdgeId<<endl;
//
//				if (intersectible(leftEdge->lower, rightEdge->lower)){
//					tmpVect.push_back(make_pair(leftEdge->EdgeId, rightEdge->EdgeId));
//					cerr<<"Vert "<<CurVert<<" Edge pair "<<leftEdge->EdgeId<<" "<<rightEdge->EdgeId<<endl;
//				}
//			}
//		}
//		if (tmpVect.size() > 1)
//			EdgePairs.insert(make_pair(CurVert,tmpVect));
//	}
//	cerr<<"Start Spliting"<<endl;
////	SplitVertecesByEdgeConnections(graph, EdgePairs, false);
//	doSplit(graph, EdgePairs);
//	cerr<<"End Spliting"<<endl;
//}
//

