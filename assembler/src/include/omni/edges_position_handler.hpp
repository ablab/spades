//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * edges_position_handler.hpp
 *
 *  Created on: 22.07.2011
 *
 */

#ifndef EDGES_POSITION_HANDLER_HPP_
#define EDGES_POSITION_HANDLER_HPP_

//#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
#include "omni_utils.hpp"
using namespace omnigraph;

namespace omnigraph {

class EdgePosition {
public:
	MappingRange m_range_;
	int start_;

	int end_;
	int start() const{return m_range_.initial_range.start_pos;}
	int end() const{return m_range_.initial_range.end_pos;}
	std::string contigId_;
	EdgePosition (int start, int end, std::string contigId = "0"):  m_range_(Range(start,end), Range(0,0)), start_(start), end_(end), contigId_(contigId) {
	};
};

bool PosCompare(const EdgePosition &a, const EdgePosition &b){
	int aend =  a.end();
	int bend =  b.end();
	return ((a.contigId_ < b.contigId_) ||((a.contigId_ == b.contigId_)&&( aend < bend )));
}

vector<EdgePosition> GluePositionsLists(vector<EdgePosition> v1, vector<EdgePosition> v2, int max_single_gap = 0){
	vector<EdgePosition> res;
	if (v1.size() == 0 && v2.size() == 0) return res;
//	if (v1.size() == 0) {res = v2;}
//	if (v2.size() == 0) {res = v1;}

	if (v1.size() == 0 || v2.size() == 0) {
		DEBUG("GluePosition fist parameter size "<<v1.size()<<", second parameter size "<<v2.size());
		set<std::string> contigs_num;
		for (auto iter = res.begin(); iter != res.end(); ++iter){
			if (contigs_num.find(iter->contigId_)==contigs_num.end()){
				DEBUG("Contig "<<iter->contigId_<< " glued with empty edge");
				contigs_num.insert(iter->contigId_);
			}
		}
		//INFO("Possible misassemble happened");
		return res;
	}
	for( size_t i = 0; i< v1.size(); i++){
		int best_fit_j = -1;
		for( size_t j = 0; j< v2.size(); j++){
			if (v1[i].contigId_ == v2[j].contigId_){
				if (v1[i].end() + 1 == v2[j].start()) {
					best_fit_j = j;
					break;
				}
				else
				{
					if ((v1[i].end() < v2[j].start())&&(v1[i].end() + max_single_gap + 1 >= v2[j].start())){
						//res.push_back(EdgePosition(v1[i].start_, v2[j].end_, v1[i].contigId_));
						if (best_fit_j < 0) best_fit_j = j;
						else if (v2[j].start() < v2[best_fit_j].start()) best_fit_j = j;
					}
				}

			}
	 	}
		if (best_fit_j != -1) {
			res.push_back(EdgePosition(v1[i].start(), v2[best_fit_j].end(), v1[i].contigId_));
			if (v2[best_fit_j].start() - v1[i].end() > 1){
				DEBUG("Contig "<<v1[i].contigId_<< " Glue parts with gap: "<<v1[i].start()<<"-"<<v1[i].end()<<" and "<<v2[best_fit_j].start()<<"-"<<v2[best_fit_j].end());
			}
		} else {
			//INFO("Possible misassemble happened");
		}
	}
	return res;
}


template<class Graph>
class EdgesPositionHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;

	int max_single_gap_;
	std::map<EdgeId, vector<EdgePosition>> EdgesPositions;
public:
	const std::map<EdgeId, vector<EdgePosition>> &edges_positions() const {
		return EdgesPositions;
	}

	const vector<EdgePosition> &GetEdgePositions(EdgeId edge) const {
		auto it = EdgesPositions.find(edge);
		VERIFY(it != EdgesPositions.end());
		return it->second;
	}

	void AddEdgePosition (EdgeId NewEdgeId, int start, int end, std::string contigId = "0") {
		if (EdgesPositions.find(NewEdgeId) == EdgesPositions.end()) {
			vector<EdgePosition> NewVec;
			EdgesPositions[NewEdgeId] = NewVec;
		}
		EdgePosition NewPos(start, end, contigId);
		(EdgesPositions[NewEdgeId]).push_back(NewPos);
//		DEBUG("Add pos "<<NewPos.start_<<" "<<NewPos.end_<<" for edge "<<NewEdgeId<<" total positions: "<< EdgesPositions[NewEdgeId].size());

		if (EdgesPositions[NewEdgeId].size()>1){
			std::sort(EdgesPositions[NewEdgeId].begin(), EdgesPositions[NewEdgeId].end(), PosCompare);
		}

	}

	void AddEdgePosition (EdgeId NewEdgeId, EdgePosition NewPos) {
		if (EdgesPositions.find(NewEdgeId) == EdgesPositions.end()) {
			vector<EdgePosition> NewVec;
			EdgesPositions[NewEdgeId] = NewVec;
		}
		(EdgesPositions[NewEdgeId]).push_back(NewPos);
	//	DEBUG("Add pos "<<NewPos.start_<<" "<<NewPos.end_<<" for edge "<<NewEdgeId<<" total positions: "<< EdgesPositions[NewEdgeId].size());

		if (EdgesPositions[NewEdgeId].size()>1){
			std::sort(EdgesPositions[NewEdgeId].begin(), EdgesPositions[NewEdgeId].end(), PosCompare);
		}

	}

	bool IsConsistentWithGenome(vector<EdgeId> Path){
		if (Path.size() > 0) {
	 		 vector<EdgePosition> res = (EdgesPositions[Path[0]]);
	 		 for (size_t i = 1; i<Path.size(); i++){
	 			 res = GluePositionsLists(res, EdgesPositions[Path[i]], max_single_gap_);
	 		 }
			if (res.size()>0){
				for(size_t i = 0; i<Path.size(); i++) {
					//todo what was it???
//					if (res[i].contigId_ < 15)
						return true; //ToDo: Curent pipeline trace genome as contigsId 0, 1, 10 and 11 but in future it can be not true.
				}
			}
		}
		return false;
	}

	void AddEdgePosition (EdgeId NewEdgeId, vector<EdgePosition> NewPositions) {
			if (EdgesPositions.find(NewEdgeId) == EdgesPositions.end()) {
				vector<EdgePosition> NewVec;
				EdgesPositions[NewEdgeId] = NewVec;
			}
			for (auto iter = NewPositions.begin(); iter != NewPositions.end(); ++iter)
				(EdgesPositions[NewEdgeId]).push_back(*iter);
		//	DEBUG("Add pos "<<NewPos.start_<<" "<<NewPos.end_<<" for edge "<<NewEdgeId<<" total positions: "<< EdgesPositions[NewEdgeId].size());

			if (EdgesPositions[NewEdgeId].size()>1){
				std::sort(EdgesPositions[NewEdgeId].begin(), EdgesPositions[NewEdgeId].end(), PosCompare);
			}

	}

	std::string str(EdgeId edgeId) const {
		std::stringstream ss;
		auto it = EdgesPositions.find(edgeId);
		if (it != EdgesPositions.end()) {
			TRACE("Number of labels " << it->second.size());
			for (auto pos_it = it->second.begin(), end = it->second.end(); pos_it != end; ++pos_it) {
				ss << "(" << pos_it->contigId_ << ": " << pos_it->start_ << "-" << pos_it->end_ << ")\\n";
			}
		}
		return ss.str();
	}

	EdgesPositionHandler(Graph &g, int max_single_gap = 0) :
		GraphActionHandler<Graph> (g, "EdgePositionHandler"), max_single_gap_(max_single_gap) {
	}

	virtual ~EdgesPositionHandler() {
		TRACE("~EdgePositionHandler ok");
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//		DEBUG("Handle glue ");

		AddEdgePosition(new_edge, (EdgesPositions[edge1]));
		AddEdgePosition(new_edge, (EdgesPositions[edge2]));

		if (EdgesPositions[edge1].size() > 0 && EdgesPositions[edge2].size() > 0) {
			DEBUG("Gluing two edges with not empty positions:");
			DEBUG("First: "<<str(edge1));
			DEBUG("Second: "<<str(edge2));
		}

/*		 for( size_t i = 0; i< EdgesPositions[edge1].size(); i++){
			 for( size_t j = 0; j< EdgesPositions[edge2].size(); j++){
//				 DEBUG(" "<<EdgesPositions[edge1])[i].start_<<" "<<EdgesPositions[edge1])[i].end_);
//				 DEBUG(" "<<EdgesPositions[edge2])[j].start_<<" "<<EdgesPositions[edge2])[j].end_);
				 if ((EdgesPositions[edge1])[i].end_ + 1 == (EdgesPositions[edge2])[j].start_) {
					 AddEdgePosition(new_edge, (EdgesPositions[edge1])[i].start_, (EdgesPositions[edge2])[j].end_);
				 }
			 }
		 }
		 */
	 }


	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		if (EdgesPositions.find(oldEdge) != EdgesPositions.end()) {
			if (EdgesPositions[oldEdge].size() > 0){
				size_t length1 = this->g().length(newEdge1);
				size_t length2 = this->g().length(newEdge2);
				for (auto iter = EdgesPositions[oldEdge].begin(); iter != EdgesPositions[oldEdge].end(); ++iter){
					int end1 = iter->start_ + (length1*(iter->end_ - iter->start_))/(length1+length2);
					AddEdgePosition(newEdge1, iter->start_, end1, iter->contigId_);
					if(iter->end_ >= end1 + 1)
						AddEdgePosition(newEdge2, end1 + 1, iter->end_, iter->contigId_);
					DEBUG("Contig "<<iter->contigId_<<" Split: " << iter->start_<<"--"<<iter->end_<<" after pos "<<end1);
				}
//				 DEBUG("EdgesPositionHandler not handled Split yet");
			}
		}
	 }

 	 virtual void HandleMerge(const vector<EdgeId>& oldEdges, EdgeId newEdge) {
//		 DEBUG("HandleMerge by position handler");
 		 // we assume that all edge have good ordered position labels.
 		 size_t n = oldEdges.size();
 		 vector<EdgePosition> res = (EdgesPositions[oldEdges[0]]);
 		 bool positive_size = (res.size() > 0);
 		 for (size_t i = 1; i < n; i++) {
 			 res = GluePositionsLists(res, EdgesPositions[oldEdges[i]], max_single_gap_);
 	 		 positive_size = positive_size || (EdgesPositions[oldEdges[i]].size() > 0);
 		 }

 		 if (positive_size && (res.size() == 0)){
 			 DEBUG("Merge operation broke some positions:");
 	 		 for (size_t i = 0; i < n; i++) {
 	 			 DEBUG("Size "<<EdgesPositions[oldEdges[i]].size()<<" "<<str(oldEdges[i]));
 	 	 		 positive_size = positive_size || (EdgesPositions[oldEdges[i]].size() > 0);
 	 		 }
 		 }

 		 AddEdgePosition(newEdge, res);
	 }
/*
	virtual void HandleAdd(VertexId v) {
		AddVertexIntId(v);
	}
	virtual void HandleDelete(VertexId v) {
		ClearVertexId(v);
	}
*/
 	virtual void HandleAdd(EdgeId e) {
 //		TRACE("Add edge "<<e);
		if (EdgesPositions.find(e) == EdgesPositions.end()) {
 			vector<EdgePosition> NewVec;
 			EdgesPositions[e] = NewVec;
		}
 	}

 	virtual void HandleDelete(EdgeId e) {
//		if (EdgesPositions[e].size() > 0) {
//			DEBUG("Delete edge "<<e<<" handled. Not empty info: "<<EdgesPositions[e].size());
//			for (size_t i = 0; i < EdgesPositions[e].size(); i++){
//				DEBUG("Position info: "<<EdgesPositions[e][i].start_<<" --- "<<EdgesPositions[e][i].end_);
//			}
//		}
//		else {
//			DEBUG("Delete edge "<<e<<" handled.");
//		}
		EdgesPositions.erase(e);
	}

	void HandleVertexSplit(VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges, vector<double> &split_coefficients, VertexId oldVertex) {
 		for (auto cur_edges_pair = newEdges.begin(); cur_edges_pair != newEdges.end(); ++cur_edges_pair){
 			AddEdgePosition(cur_edges_pair->second, EdgesPositions[cur_edges_pair->first]);
 		}
 	}

private:
    DECL_LOGGER("EdgesPositionHandler");
};

}

#endif /* EDGES_POSITION_HANDLER_HPP_ */
