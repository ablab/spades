//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "abstract_editable_graph.hpp"
#include "id_track_handler.hpp"

namespace omnigraph {


class CoveredEdge {
private:
	int coverage_;

public:
	CoveredEdge() :
			coverage_(0) { }

	void SetCoverage(int coveradge) {
		coverage_ = coveradge;
	}

	void IncCoverage(int value) {
		coverage_ += value;
	}

	//not length normalized
	int GetRawCoverage() const {
		return coverage_;
	}
};



template<typename VertexIdT, typename EdgeIdT, class DataMasterT>
class AbstractGraph: public AbstractEditableGraph<VertexIdT, EdgeIdT, DataMasterT, typename set<VertexIdT>::const_iterator> {
	typedef AbstractEditableGraph<VertexIdT, EdgeIdT, DataMasterT, typename set<VertexIdT>::const_iterator> base;

public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;
	typedef DataMasterT DataMaster;
	typedef typename DataMaster::VertexData VertexData;
	typedef typename DataMaster::EdgeData EdgeData;
	typedef typename base::VertexIterator VertexIterator;


protected:

	set<VertexId> vertices_;

public:
	AbstractGraph(HandlerApplier<VertexId, EdgeId>* applier,
			const DataMaster& master) :
		base(applier, master) {
	}

	virtual ~AbstractGraph() {
		TRACE("~AbstractGraph");
	}

	virtual const vector<EdgeId> OutgoingEdges(VertexId v) const {
		return v->OutgoingEdges();
	}

	virtual const vector<EdgeId> IncomingEdges(VertexId v) const {
		return v->IncomingEdges();
	}

	virtual size_t OutgoingEdgeCount(VertexId v) const {
		return v->OutgoingEdgeCount();
	}

	virtual size_t IncomingEdgeCount(VertexId v) const {
		return v->IncomingEdgeCount();
	}

	virtual vector<EdgeId> GetEdgesBetween(VertexId v, VertexId u) const {
		return v->OutgoingEdgesTo(u);
	}

	virtual const EdgeData& data(EdgeId edge) const {
		return edge->data();
	}

	virtual const VertexData& data(VertexId v) const {
		return v->data();
	}

	virtual VertexIterator begin() const {
		return vertices_.begin();
	}

	virtual VertexIterator end() const {
		return vertices_.end();
	}

	size_t size() const {
		return vertices_.size();
	}

	virtual VertexId EdgeStart(EdgeId edge) const {
		return edge->start();
	}

	virtual VertexId EdgeEnd(EdgeId edge) const {
		return edge->end();
	}


protected:

	void AddVertexToGraph(VertexId vertex) {
		auto result = vertices_.insert(vertex);
		VERIFY(result.second); // was not in set before
	}

	void DeleteVertexFromGraph(VertexId vertex) {
		auto it = vertices_.find(vertex);
		VERIFY(it != vertices_.end()); // is it in set
		vertices_.erase(it);
	}

private:
DECL_LOGGER("AbstractGraph");
};
}
