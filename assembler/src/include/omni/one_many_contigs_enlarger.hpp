//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * one-many_contigs_enlarger.hpp
 *
 *  Created on: Aug 1, 2011
 *      Author: undead
 */

#ifndef ONE_MANY_CONTIGS_ENLARGER_HPP_
#define ONE_MANY_CONTIGS_ENLARGER_HPP_

using namespace omnigraph;

template <class Graph>
class one_many_contigs_enlarger {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	Graph &g_;
	int insert_size_;
public:
	one_many_contigs_enlarger(Graph &g, int insert_size): g_(g), insert_size_(insert_size){
	}
	void Loops_resolve(){
		return; //temporary not used
//ToDo: Think about good loop resolver
		INFO("----Resolving Loops----");
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter){
			if ((g_.EdgeStart(*iter) == g_.EdgeEnd(*iter))&&(g_.length(*iter) < insert_size_ * 1.1)){
				VertexId vertex = g_.EdgeEnd(*iter);
				EdgeId loopEdge = *iter;
				if ((g_.OutgoingEdgeCount(vertex) == 2)&&(g_.IncomingEdgeCount(vertex) == 2)){
					VertexId newVertex = g_.AddVertex();
					vector<EdgeId> outEdges = g_.OutgoingEdges(vertex);
					vector<EdgeId> inEdges = g_.IncomingEdges(vertex);
					EdgeId outEdge = (outEdges[0] == *iter)? outEdges[1]:outEdges[0];
					EdgeId inEdge = (inEdges[0] == *iter)? inEdges[1]:inEdges[0];
					EdgeId newOutEdge = g_.AddEdge(newVertex, g_.EdgeEnd(outEdge), g_.EdgeNucls(outEdge));
					EdgeId newLoopEdge = g_.AddEdge(vertex, newVertex, g_.EdgeNucls(loopEdge));
					vector<pair<EdgeId, EdgeId>> cloneEdges= {make_pair(outEdge, newOutEdge), make_pair(loopEdge, newLoopEdge)};
					vector<double> coeff = {1.,1.};
					g_.FireVertexSplit(newVertex, cloneEdges, coeff, vertex);
					g_.DeleteEdge(loopEdge);
					g_.DeleteEdge(outEdge);
					vector<EdgeId> toMerge = {inEdge, newLoopEdge, newOutEdge};
					g_.MergePath(toMerge);
				}
				else {
					g_.DeleteEdge(*iter);
				}
			}
		}
	}
/*
	void one_many_resolve(){
		INFO("one_many_resolve");
		Loops_resolve();
		int inc_count;
		int out_count;
		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter) {
			VertexId vertex = *iter;
			DEBUG(vertex);
			if ((g_.OutgoingEdgeCount(vertex) == 1) && ((inc_count = g_.IncomingEdgeCount(vertex)) >= 1)){
				EdgeId unique = g_.GetUniqueOutgoingEdge(vertex);
				if ((g_.EdgeStart(unique) != g_.EdgeEnd(unique))){
					vector<EdgeId> incEdges = g_.IncomingEdges(vertex);
					for(int j = 0; j < inc_count; j++) {
						VertexId tmp_v = g_.AddVertex();
						EdgeId edge2 = g_.AddEdge(tmp_v, g_.EdgeEnd(unique), g_.EdgeNucls(unique));
						EdgeId edge1 = g_.AddEdge(g_.EdgeStart(incEdges[j]), tmp_v, g_.EdgeNucls(incEdges[j]));
						vector<EdgeId> toMerge;
						toMerge.push_back(edge1);
						toMerge.push_back(edge2);
						DEBUG("first part ");
						g_.MergePath(toMerge);
					}
					g_.ForceDeleteVertex(vertex);
				}
			}
		}


		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter){
			VertexId vertex = *iter;
			if (((out_count = g_.OutgoingEdgeCount(vertex)) >= 1) && (g_.IncomingEdgeCount(vertex) == 1)){
				EdgeId unique = g_.GetUniqueIncomingEdge(vertex);
				if ((g_.EdgeStart(unique) != g_.EdgeEnd(unique))) {
					vector<EdgeId> outEdges = g_.OutgoingEdges(vertex);
					for(int j = 0; j < out_count; j++) {
						VertexId tmp_v = g_.AddVertex();
						EdgeId edge2 = g_.AddEdge(g_.EdgeStart(unique), tmp_v, g_.EdgeNucls(unique));
						EdgeId edge1 = g_.AddEdge(tmp_v, g_.EdgeEnd(outEdges[j]), g_.EdgeNucls(outEdges[j]));
						vector<EdgeId> toMerge ;
						toMerge.push_back(edge2);
						toMerge.push_back(edge1);
						//					DEBUG("second part ");
						g_.MergePath(toMerge);
					}
					g_.ForceDeleteVertex(vertex);
					//				DEBUG("second vertex deleted ");
				}
			}
		}



	}

*/
	void one_many_resolve_with_vertex_split(){
		INFO("one_many_resolve");
//		Loops_resolve();
		int inc_count;
		int out_count;

		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter) {
			VertexId vertex = *iter;
			DEBUG(vertex);
			if ((g_.OutgoingEdgeCount(vertex) == 1) && ((inc_count = g_.IncomingEdgeCount(vertex)) >= 1)){
				EdgeId unique = g_.GetUniqueOutgoingEdge(vertex);
				if ((g_.EdgeStart(unique) != g_.EdgeEnd(unique))){
					vector<EdgeId> incEdges = g_.IncomingEdges(vertex);
					for(int j = 0; j < inc_count; j++) {
						vector<EdgeId> SplitVect;
						SplitVect.push_back(unique);
						SplitVect.push_back(incEdges[j]);
						pair<VertexId, vector<pair<EdgeId, EdgeId>>> tmp_pair = g_.SplitVertex(vertex, SplitVect);
						VertexId tmp_v = tmp_pair.first;
						EdgeId edge2 = g_.GetUniqueOutgoingEdge(tmp_v);
						EdgeId edge1 = g_.GetUniqueIncomingEdge(tmp_v);
						vector<EdgeId> toMerge;
						toMerge.push_back(edge1);
						toMerge.push_back(edge2);

						DEBUG("first part ");
						g_.MergePath(toMerge);

					}

					g_.ForceDeleteVertex(vertex);
				}
			}
		}


		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter){
			VertexId vertex = *iter;
			if (((out_count = g_.OutgoingEdgeCount(vertex)) >= 1) && (g_.IncomingEdgeCount(vertex) == 1)){
				EdgeId unique = g_.GetUniqueIncomingEdge(vertex);
				if ((g_.EdgeStart(unique) != g_.EdgeEnd(unique))) {
					vector<EdgeId> outEdges = g_.OutgoingEdges(vertex);
					for(int j = 0; j < out_count; j++) {
						vector<EdgeId> SplitVect;
						SplitVect.push_back(unique);
						SplitVect.push_back(outEdges[j]);
						pair<VertexId, vector<pair<EdgeId, EdgeId>>> tmp_pair = g_.SplitVertex(vertex, SplitVect);
						VertexId tmp_v = tmp_pair.first;
						EdgeId edge2 = g_.GetUniqueOutgoingEdge(tmp_v);
						EdgeId edge1 = g_.GetUniqueIncomingEdge(tmp_v);
						vector<EdgeId> toMerge;
						toMerge.push_back(edge1);
						toMerge.push_back(edge2);
						DEBUG("first part ");
						g_.MergePath(toMerge);

					}
					g_.ForceDeleteVertex(vertex);
				}

			}
		}
	}
};

#endif /* ONE_MANY_CONTIGS_ENLARGER_HPP_ */
