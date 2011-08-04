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
public:
	one_many_contigs_enlarger(Graph &g): g_(g){
	}
	void one_many_resolve(){
		INFO("one_many_resolve");
		int inc_count;
		int out_count;
		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter) {
			VertexId vertex = *iter;
			DEBUG(vertex);
			if ((g_.OutgoingEdgeCount(vertex) == 1) && ((inc_count = g_.IncomingEdgeCount(vertex)) > 1)){
				EdgeId unique = g_.GetUniqueOutgoingEdge(vertex);
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


		for (auto iter = g_.SmartVertexBegin(); ! iter.IsEnd(); ++iter){
			VertexId vertex = *iter;
			if (((out_count = g_.OutgoingEdgeCount(vertex)) > 1) && (g_.IncomingEdgeCount(vertex) == 1)){
				EdgeId unique = g_.GetUniqueIncomingEdge(vertex);
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
};

#endif /* ONE_MANY_CONTIGS_ENLARGER_HPP_ */
