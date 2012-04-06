#pragma once

#include "abstract_editable_graph.hpp"
#include "id_track_handler.hpp"

namespace omnigraph {
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

private:
	typedef SmartSet<AbstractGraph, VertexId> Vertices;

	Vertices smart_vertices_;

protected:

	const set<VertexId> &vertices_;

public:
	AbstractGraph(HandlerApplier<VertexId, EdgeId>* applier,
			const DataMaster& master) :
		base(applier, master), smart_vertices_(*this), vertices_(smart_vertices_.inner_set()) {
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

private:
DECL_LOGGER("AbstractGraph");
};
}
