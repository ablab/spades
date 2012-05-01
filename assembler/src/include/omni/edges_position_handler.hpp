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
	int m_start() const{return m_range_.mapped_range.start_pos;}
	int m_end() const{return m_range_.mapped_range.end_pos;}

	std::string contigId_;
	EdgePosition (int start, int end, std::string contigId = "0", int mapped_start = 0, int mapped_end = 0):  m_range_(Range(start,end), Range(mapped_start,mapped_end)), start_(start), end_(end), contigId_(contigId) {
	};
//	EdgePosition (MappingRange& m_range, std::string contigId = "0"):  m_range_(m_range), start_(m_range.initial_range.start_pos), end_(m_range.initial_range.end_pos), contigId_(contigId) {
//	};
	EdgePosition (MappingRange& m_range, std::string contigId = "0", int shift = 0):  m_range_(m_range), start_(m_range.initial_range.start_pos), end_(m_range.initial_range.end_pos), contigId_(contigId) {
		m_range.mapped_range.start_pos += shift;
		m_range.mapped_range.end_pos += shift;
	};
//	EdgePosition (EdgePosition& e_pos, int shift = 0):  m_range_(e_pos.m_range_), start_(m_range_.initial_range.start_pos), end_(m_range_.initial_range.end_pos), contigId_(e_pos.contigId_) {
//		m_range_.mapped_range.start_pos += shift;
//		m_range_.mapped_range.end_pos += shift;
//	};

	EdgePosition (const EdgePosition& e_pos, int shift = 0):  m_range_(e_pos.m_range_), start_(m_range_.initial_range.start_pos), end_(m_range_.initial_range.end_pos), contigId_(e_pos.contigId_) {
		m_range_.mapped_range.start_pos += shift;
		m_range_.mapped_range.end_pos += shift;
	};

	void shift_mapped_range(int shift){
		if (m_range_.mapped_range.end_pos != 0) {
			m_range_.mapped_range.start_pos += shift;
			m_range_.mapped_range.end_pos += shift;
		}
	}
};

bool PosCompare(const EdgePosition &a, const EdgePosition &b){
	int aend =  a.end();
	int bend =  b.end();
	return ((a.contigId_ < b.contigId_) ||((a.contigId_ == b.contigId_)&&( aend < bend )));
}

vector<EdgePosition> RangeGluePositionsLists(vector<EdgePosition> v1, vector<EdgePosition> v2, int max_single_gap = 0, int shift = 0){
	DEBUG("RangeGluePositionsLists: v1.size = "<<v1.size()<<" v2.size = "<< v2.size());
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
		if (v2.size() == 0) {res = v1;}
		if (v1.size() == 0) {
			for (size_t i = 0; i < v2.size(); i++){
				res.push_back(EdgePosition(v2[i], shift));
			}
		}

		return res;
	}
	set<size_t> used_seconds;
	for( size_t i = 0; i< v1.size(); i++){
		int best_fit_j = -1;
		for( size_t j = 0; j< v2.size(); j++){
			if (v1[i].contigId_ == v2[j].contigId_){
				if ((v1[i].end() + 1 == v2[j].start()) && (v1[i].m_end() + 1 == v2[j].m_start() + shift)) {
					best_fit_j = j;
					break;
				}
				else
				{
					if ((v1[i].end() < v2[j].start())&&(v1[i].end() + max_single_gap + 1 >= v2[j].start())
					 && (v1[i].m_end() + max_single_gap + 1 >= v2[j].m_start() + shift)){
						//res.push_back(EdgePosition(v1[i].start_, v2[j].end_, v1[i].contigId_));
						if (best_fit_j < 0) best_fit_j = j;
						else if (v2[j].start() < v2[best_fit_j].start()) best_fit_j = j;
					}
				}

			}
	 	}
		if (best_fit_j != -1) {
			res.push_back(EdgePosition(v1[i].start(), v2[best_fit_j].end(), v1[i].contigId_, v1[i].m_start(), v2[best_fit_j].m_end() + shift));
			used_seconds.insert(best_fit_j);
			if (v2[best_fit_j].start() - v1[i].end() > 1){
				DEBUG("Contig "<<v1[i].contigId_<< " Glue parts with gap: "<<v1[i].start()<<"-"<<v1[i].end()<<" and "<<v2[best_fit_j].start()<<"-"<<v2[best_fit_j].end());
			}
		} else {
			//INFO("Possible misassemble happened");
			res.push_back(v1[i]);
		}
	}
	for( size_t j = 0; j< v2.size(); j++){
		if (used_seconds.find(j) == used_seconds.end()){
			res.push_back(EdgePosition(v2[j], shift));
		}
	}
	return res;
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
	bool careful_ranges_;
public:
	bool is_careful(){return careful_ranges_;};
	const std::map<EdgeId, vector<EdgePosition>> &edges_positions() const {
		return EdgesPositions;
	}

	const vector<EdgePosition> &GetEdgePositions(EdgeId edge) const {
		auto it = EdgesPositions.find(edge);
		VERIFY(it != EdgesPositions.end());
		return it->second;
	}

	void AddEdgePosition (EdgeId NewEdgeId, int start, int end, std::string contigId = "0", int m_start = 0, int m_end = 0) {
		if (EdgesPositions.find(NewEdgeId) == EdgesPositions.end()) {
			vector<EdgePosition> NewVec;
			EdgesPositions[NewEdgeId] = NewVec;
		}
		if (!careful_ranges_){
			m_start = 0;
			m_end = 0;
		}
		EdgePosition NewPos(start, end, contigId, m_start, m_end);
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
				if (careful_ranges_){
					ss << "(" << pos_it->contigId_ << ": " << pos_it->start() << "-" << pos_it->end() << ": "<< pos_it->m_start() << "-" << pos_it->m_end() << ")\\n";
				}
				else {
					ss << "(" << pos_it->contigId_ << ": " << pos_it->start_ << "-" << pos_it->end_ << ")\\n";
				}

			}
		}
		return ss.str();
	}

	EdgesPositionHandler(Graph &g, int max_single_gap = 0, bool careful_ranges = false) :
		GraphActionHandler<Graph> (g, "EdgePositionHandler"), max_single_gap_(max_single_gap), careful_ranges_(careful_ranges) {
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
				if (careful_ranges_){
					int length1 = this->g().length(newEdge1);
//					int length2 = this->g().length(newEdge2);
					for (auto iter = EdgesPositions[oldEdge].begin(); iter != EdgesPositions[oldEdge].end(); ++iter){
						if (iter->m_end() <= length1){
							AddEdgePosition(newEdge1, iter->start(), iter->end(), iter->contigId_, iter->m_start(), iter->m_end());
						} else
						if (iter->m_start() > length1){
							AddEdgePosition(newEdge2, iter->start(), iter->end(), iter->contigId_, iter->m_start() - length1, iter->m_end() - length1);
						} else {
							int segm1len =  (iter->end() - iter->start() + 1) * (length1 - iter->m_start() + 1) / (iter->m_end() - iter->m_start() + 1);
							int segm2len =  iter->m_end() - iter->m_start() - segm1len + 1;
							if (segm1len > 0)
								AddEdgePosition(newEdge1, iter->start(), iter->start() + segm1len - 1, iter->contigId_, iter->m_start(), length1);
							if (segm2len > 0)
								AddEdgePosition(newEdge2, iter->start() + segm1len, iter->end(), iter->contigId_, 1, iter->m_end() - length1);
						}
//						int end1 = iter->start_ + (length1*(iter->end_ - iter->start_))/(length1+length2);
//						AddEdgePosition(newEdge1, iter->start_, end1, iter->contigId_);
//						if(iter->end_ >= end1 + 1)
//							AddEdgePosition(newEdge2, end1 + 1, iter->end_, iter->contigId_);
//						DEBUG("Contig "<<iter->contigId_<<" Split: " << iter->start_<<"--"<<iter->end_<<" after pos "<<end1);
					}
				}
				else {
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
	 }

 	 virtual void HandleMerge(const vector<EdgeId>& oldEdges, EdgeId newEdge) {
//		 DEBUG("HandleMerge by position handler");
 		 // we assume that all edge have good ordered position labels.
 		 size_t n = oldEdges.size();
 		 vector<EdgePosition> res = (EdgesPositions[oldEdges[0]]);
 		 int shift = this->g().length(oldEdges[0]);
 		 bool positive_size = (res.size() > 0);
 		 for (size_t i = 1; i < n; i++) {
 			 if (careful_ranges_){
 				 res = RangeGluePositionsLists(res, EdgesPositions[oldEdges[i]], max_single_gap_, shift);
 		 		 shift += this->g().length(oldEdges[i]);
 			 }
 			 else {
 				 res = GluePositionsLists(res, EdgesPositions[oldEdges[i]], max_single_gap_);
 			 }

 			 positive_size = positive_size || (EdgesPositions[oldEdges[i]].size() > 0);
 		 }

 		 if (positive_size && (res.size() == 0)){
 			 DEBUG("Merge operation broke some positions:");
 	 		 for (size_t i = 0; i < n; i++) {
 	 			 DEBUG("Size "<<EdgesPositions[oldEdges[i]].size()<<" "<<str(oldEdges[i]));
 	 	 		 positive_size = positive_size || (EdgesPositions[oldEdges[i]].size() > 0);
 	 		 }
 		 }
 		 DEBUG("NEW res.size = "<< res.size());
 		 AddEdgePosition(newEdge, res);
	 }

 	 virtual void HandleAdd(EdgeId e) {
 //		TRACE("Add edge "<<e);
		if (EdgesPositions.find(e) == EdgesPositions.end()) {
 			vector<EdgePosition> NewVec;
 			EdgesPositions[e] = NewVec;
		}
 	}

 	virtual void HandleDelete(EdgeId e) {
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
