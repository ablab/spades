/*
 * total_labeler.hpp
 *
 *  Created on: 10.08.2011
 */

#ifndef TOTAL_LABELER_HPP_
#define TOTAL_LABELER_HPP_

#include "graph_labeler.hpp"
#include "simple_tools.hpp"
#include "edge_labels_handler.hpp"
#include "ID_track_handler.hpp"
#include "edges_position_handler.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <map>

using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class TotalLabelerGraphStruct {
public:
	Graph& g_;
	IdTrackHandler<Graph>* IDs;
	EdgesPositionHandler<Graph>* EdgesPos;
	EdgeLabelHandler<Graph>* EdgesLabels;
	TotalLabelerGraphStruct(Graph &g, IdTrackHandler<Graph>* id_handler, EdgesPositionHandler<Graph>* position_handler,
							EdgeLabelHandler<Graph>* label_handler = NULL): g_(g),
							IDs(id_handler), EdgesPos(position_handler), EdgesLabels(label_handler)	{

	}
	~TotalLabelerGraphStruct(){};
};

template<class Graph>
class TotalLabeler: public GraphLabeler<Graph> {

protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	TotalLabelerGraphStruct<Graph>* graph_struct;
	TotalLabelerGraphStruct<Graph>* proto_graph_struct;

public:
	TotalLabeler(TotalLabelerGraphStruct<Graph>* g_struct, TotalLabelerGraphStruct<Graph>* proto_g_struct = NULL) :
		graph_struct(g_struct), proto_graph_struct(proto_g_struct)  {
	}

	virtual std::string label(VertexId vertexId) const {
		return graph_struct->g_.str(vertexId);
	}

	std::string OldEdgeIdToStr(EdgeId e_id){
		if (proto_graph_struct != NULL)
			if (proto_graph_struct->IDs != NULL){
				int id = proto_graph_struct->IDs->ReturnIntId(e_id);
				return ToString(id);
			}
		return ToString(e_id);
	}

	virtual std::string label(EdgeId edgeId) const {
		std::string ret_label;
		if (graph_struct->IDs != NULL) {
			ret_label += "Id "+graph_struct->IDs->str(edgeId)+"\\n";
		}
		if (graph_struct->EdgesPos != NULL){
			ret_label += "Positions:\\n"+ graph_struct->EdgesPos->str(edgeId);
		}
		if (graph_struct->EdgesLabels != NULL){
			if ((proto_graph_struct != NULL) && (proto_graph_struct->IDs != NULL)) {
				boost::function<string (EdgeId)> f = boost::bind(&IdTrackHandler<Graph>::str, *(proto_graph_struct->IDs), _1);
				ret_label += "Labels:\\n" + graph_struct->EdgesLabels->str(edgeId, f);
			}
			else {
				ret_label += "Labels:\\n" + graph_struct->EdgesLabels->str(edgeId);
			}
		}
		ret_label += "Len(cov): " + graph_struct->g_.str(edgeId);
		return ret_label;
	}
	virtual ~TotalLabeler() {
		TRACE("~TotalLabeler");
	}

};

}


#endif /* TOTAL_LABELER_HPP_ */
