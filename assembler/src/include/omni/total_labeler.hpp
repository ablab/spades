//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
#include "id_track_handler.hpp"
#include "edges_position_handler.hpp"

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <map>

using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class TotalLabelerGraphStruct
    : boost::noncopyable
{
public:
	const Graph& g_;
	const IdTrackHandler<Graph>* IDs;
	const EdgesPositionHandler<Graph>* EdgesPos;
	const EdgeLabelHandler<Graph>* EdgesLabels;
	TotalLabelerGraphStruct(const Graph &g, const IdTrackHandler<Graph>* id_handler, const EdgesPositionHandler<Graph>* position_handler,
			const EdgeLabelHandler<Graph>* label_handler = NULL): g_(g),
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
public:
	const TotalLabelerGraphStruct<Graph>* graph_struct;
	const TotalLabelerGraphStruct<Graph>* proto_graph_struct;

	TotalLabeler(const TotalLabelerGraphStruct<Graph>* g_struct, const TotalLabelerGraphStruct<Graph>* proto_g_struct = NULL) :
		graph_struct(g_struct), proto_graph_struct(proto_g_struct)  {

		if ((proto_graph_struct != NULL)){
			VERIFY(proto_graph_struct->g_.ReturnIntIdPointer());
		}
	}

	virtual std::string label(VertexId vertexId) const {
		int vId = graph_struct->IDs->ReturnIntId(vertexId);
		return ToString(vId);
	}

	std::string OldEdgeIdToStr(EdgeId e_id) const {
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
				boost::function<string (EdgeId)> f = boost::bind(&IdTrackHandler<Graph>::str, boost::ref(*(proto_graph_struct->IDs)), _1);
				ret_label += "Labels:\\n" + graph_struct->EdgesLabels->str(edgeId, f);
			}
			else {
				ret_label += "Labels:\\n" + graph_struct->EdgesLabels->str(edgeId);
			}
		}



		int len = graph_struct->g_.length(edgeId);

		double cov = graph_struct->g_.coverage(edgeId);

		ret_label += "Len(cov): " + ToString(len)+"("+ToString(cov)+")";// + graph_struct->g_.str(edgeId);

		return ret_label;
	}

	virtual ~TotalLabeler() {
		TRACE("~TotalLabeler");
	}

};

}


#endif /* TOTAL_LABELER_HPP_ */
