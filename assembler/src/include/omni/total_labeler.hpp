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

#include "visualization/graph_labeler.hpp"
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
class TotalLabeler: public GraphLabeler<Graph> {
private:
    const Graph& g_;
    const EdgesPositionHandler<Graph> &edges_positions_;
protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
public:

	TotalLabeler(const Graph &g, const EdgesPositionHandler<Graph> &position_handler) :
		g_(g), edges_positions_(position_handler) {
	}

	virtual std::string label(VertexId vertexId) const {
		return ToString(vertexId.int_id());
	}

	std::string edge_id_str(EdgeId e) const {
	    return graph_struct->g_.str(e);
	}

	virtual std::string label(EdgeId edgeId) const {
		std::string ret_label;
		ret_label += "Id "+graph_struct->g_.str(edgeId)+"\\n";
		ret_label += "Positions:\\n"+ graph_struct->EdgesPos->str(edgeId);
		size_t len = g_.length(edgeId);
		double cov = g_.coverage(edgeId);
		ret_label += "Len(cov): " + ToString(len) + "(" + ToString(cov) + ")";  // + graph_struct->g_.str(edgeId);
		return ret_label;
	}

	virtual ~TotalLabeler() {
		TRACE("~TotalLabeler");
	}

};

}


#endif /* TOTAL_LABELER_HPP_ */
