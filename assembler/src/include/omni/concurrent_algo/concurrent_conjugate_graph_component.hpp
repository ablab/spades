//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
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
	typedef typename Graph::VertexData VertexData;

	template<class InputVertexIterator>
	ConcurrentConjugateGraphComponent(
			Graph& graph,
			const restricted::PeriodicIdDistributor& id_distributor,
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
		return this->graph_.conjugate(vertex);
	}

	EdgeId conjugate(EdgeId edge) const {
		return this->graph_.conjugate(edge);
	}

	virtual bool IsInternalSafe(const VertexId& vertex) const {
		return this->IsInternal(vertex) && this->IsInternal(conjugate(vertex));
	}

	virtual bool IsInternalSafe(const EdgeId& edge) const {
		return
				IsInternalSafe(this->EdgeStart(edge)) &&
				IsInternalSafe(this->EdgeEnd(edge));
	}

	virtual bool IsInComponentSafe(const EdgeId& edge) const {
		return
			this->IsInComponent(edge) &&
			this->IsInComponent(this->graph_.conjugate(edge));
	}

	virtual ~ConcurrentConjugateGraphComponent() {
	}


protected:

	virtual void AddVertexToComponent(VertexId vertex) {
		this->vertices_.insert(vertex);
		this->temporary_vertices_.insert(vertex);

		this->vertices_.insert(GetConjugateWithoutChecks(vertex));
		this->temporary_vertices_.insert(GetConjugateWithoutChecks(vertex));
	}

	virtual VertexId HiddenAddVertex(const VertexData &data) {
		VertexId vertex = this->CreateVertex(data);
		AddVertexToComponent(vertex);
		return vertex;
	}

	virtual void HiddenDeleteVertex(VertexId vertex) {
//		VERIFY(IsInternalSafe(vertex));

		VertexId conjugate_vertex = conjugate(vertex);

		this->vertices_.erase(vertex);
		this->vertices_.erase(conjugate_vertex);

		if (this->temporary_vertices_.find(vertex) != this->temporary_vertices_.end()) {
			this->temporary_vertices_.erase(vertex);
			this->temporary_vertices_.erase(conjugate_vertex);

			this->DestroyVertex(vertex); // conjugate will be deleted too
		} else {
			this->deleted_vertices_.push_back(vertex);
		}
	}

	VertexId GetConjugateWithoutChecks(VertexId vertex) const {
		return this->graph_.conjugate(vertex);
	}

	EdgeId GetConjugateWithoutChecks(EdgeId edge) const {
		return this->graph_.conjugate(edge);
	}
};

} //namespace omnigraph

#endif /* CONCURRENT_CONJUGATE_GRAPH_COMPONENT_HPP_ */
