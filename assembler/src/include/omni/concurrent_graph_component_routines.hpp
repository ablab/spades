//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * concurrent_graph_component_routines.hpp
 *
 *  Created on: Aug 14, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef CONCURRENT_GRAPH_COMPONENT_ROUTINES_HPP_
#define CONCURRENT_GRAPH_COMPONENT_ROUTINES_HPP_


#include <boost/utility.hpp>

#include <exception>

#include "omni/order_and_law.hpp"
#include "omni/omni_utils.hpp"


namespace omnigraph {


class PoolEdgeIdDistributor : public restricted::IdDistributor {

public:
	PoolEdgeIdDistributor(size_t pool_size) {
		cur_id_ = restricted::GlobalIdDistributor::GetInstance()->GetIdPool(pool_size);
		max_id_ = cur_id_ + pool_size - 1;
	}

	virtual size_t GetId() {
		VERIFY(HasIds() > 0);
		return cur_id_++;
	}

	size_t HasIds() {
		return max_id_ - cur_id_;
	}

private:
	size_t cur_id_;
	size_t max_id_;
};


template<typename VertexId, typename EdgeId>
class HandlerApplierStub : public omnigraph::HandlerApplier<VertexId, EdgeId> {
public:

	virtual void
	ApplyAdding(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const {
		VERIFY(false);
	}

	virtual void
	ApplyAdding(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const {
		VERIFY(false);
	}

	virtual void
	ApplyAdd(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const {
		VERIFY(false);
	}

	virtual void
	ApplyAdd(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const {
		VERIFY(false);
	}

	virtual void
	ApplyDelete(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const {
		VERIFY(false);
	}

	virtual void
	ApplyDelete(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const {
		VERIFY(false);
	}

	virtual void ApplyMerge(ActionHandler<VertexId, EdgeId> *handler,
	vector<EdgeId> old_edges, EdgeId new_edge) const {
		VERIFY(false);
	}

	virtual void ApplyGlue(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
		VERIFY(false);
	}

	virtual void ApplySplit(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const {
		VERIFY(false);
	}

	virtual void ApplyVertexSplit(ActionHandler<VertexId, EdgeId> *handler,
	VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges,
	vector<double> &split_coefficients, VertexId oldVertex) const {
		VERIFY(false);
	}

	virtual ~HandlerApplierStub() {
	}
};


} // namespace debuijn


#endif /* CONCURRENT_GRAPH_COMPONENT_ROUTINES_HPP_ */
