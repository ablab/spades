#pragma once

#include "abstract_editable_graph.hpp"
#include "id_track_handler.hpp"

namespace omnigraph {
template<typename VertexIdT, typename EdgeIdT, class DataMasterT,
		typename VertexIt>
class AbstractGraph: public AbstractEditableGraph<VertexIdT, EdgeIdT,
		DataMasterT, VertexIt> {
	typedef AbstractEditableGraph<VertexIdT, EdgeIdT, DataMasterT, VertexIt>
			base;
public:
	typedef VertexIdT VertexId;
	typedef EdgeIdT EdgeId;
	typedef DataMasterT DataMaster;
	typedef typename DataMaster::VertexData VertexData;
	typedef typename DataMaster::EdgeData EdgeData;
	typedef VertexIt VertexIterator;

protected:
	typedef set<VertexId> Vertices;

	Vertices vertices_;

	mutable BaseIdTrackHandler<VertexIdT, EdgeIdT>* int_ids_;

public:
	AbstractGraph(HandlerApplier<VertexId, EdgeId>* applier,
			const DataMaster& master) :
		base(applier, master), int_ids_(NULL){
	}

	virtual ~AbstractGraph() {
		TRACE("~AbstractGraph");
	}

	void set_int_ids(BaseIdTrackHandler<VertexIdT, EdgeIdT>* int_ids) const {
		VERIFY(!int_ids_ || !int_ids);
		int_ids_ = int_ids;
	}

	BaseIdTrackHandler<VertexIdT, EdgeIdT> &int_ids() const {
		VERIFY(int_ids_);
		return *int_ids_;
	}

	int int_id(EdgeId edge) const {
		VERIFY(int_ids_);
		return int_ids_->ReturnIntId(edge);
	}

	int int_id(VertexId vertex) const {
		VERIFY(int_ids_);
		return int_ids_->ReturnIntId(vertex);
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
