//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************


/*
 * concurrent_conjugate_graph_component.hpp
 *
 *  Created on: Aug 20, 2012
 *      Author: Alexander Opeykin (alexander.opeykin@gmail.com)
 */


#ifndef CONCURRENT_CONJUGATE_GRAPH_COMPONENT_HPP_
#define CONCURRENT_CONJUGATE_GRAPH_COMPONENT_HPP_

#include "concurrent_graph_component.hpp"
#include "omni_utils.hpp"

namespace omnigraph {

template <typename Graph>
class ConcurrentConjugateGraphComponent : public ConcurrentGraphComponent<Graph> {

public:
	typedef ConcurrentGraphComponent<Graph> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	template<class InputVertexIterator>
	ConcurrentConjugateGraphComponent(
			Graph& graph,
			restricted::IdDistributor& id_distributor,
			InputVertexIterator verticesBegin,
			InputVertexIterator verticesEnd)
				: base(
						graph,
						new PairedHandlerApplier<ConcurrentConjugateGraphComponent>(*this),
						id_distributor,
						verticesBegin,
						verticesEnd) {
	}

	VertexId conjugate(VertexId vertex) const {
		this->SetFlagIfNotInComponent(vertex);

		VertexId conjugate = this->graph_.conjugate(vertex);
		this->SetFlagIfNotInComponent(conjugate);

		return conjugate;
	}

	EdgeId conjugate(EdgeId edge) const {
		this->SetFlagIfNotInComponent(edge);

		EdgeId conjugate = this->graph_.conjugate(edge);
		this->SetFlagIfNotInComponent(conjugate);

		return conjugate;
	}

	virtual bool IsInternalSafe(const VertexId& vertex) const {
		return this->IsInternal(vertex) && this->IsInternal(conjugate(vertex));
	}

	virtual bool IsInternalSafe(const EdgeId& edge) const {
		EdgeId conjugate = this->conjugate(edge);

		return
				IsInternalSafe(this->EdgeStart(edge)) &&
				IsInternalSafe(this->EdgeEnd(edge)) &&
				IsInternalSafe(this->EdgeStart(conjugate)) &&
				IsInternalSafe(this->EdgeEnd(conjugate));
	}

	virtual ~ConcurrentConjugateGraphComponent() {
	}

protected:
	virtual void HiddenDeleteVertex(VertexId vertex) {
		VertexId conjugate = this->graph_.conjugate(vertex);

		VERIFY(this->all_actions_valid_);
		VERIFY(this->IsInComponent(vertex));
		VERIFY(this->IsInComponent(conjugate));

		VERIFY(this->IncomingEdgeCount(vertex) == 0);
		VERIFY(this->OutgoingEdgeCount(vertex) == 0);

		VERIFY(this->IncomingEdgeCount(conjugate) == 0);
		VERIFY(this->OutgoingEdgeCount(conjugate) == 0);

		this->vertices_.erase(vertex);
		this->vertices_.erase(conjugate);
	}

};

} //namespace omnigraph

#endif /* CONCURRENT_CONJUGATE_GRAPH_COMPONENT_HPP_ */
